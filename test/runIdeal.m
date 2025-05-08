% if you acquired this file from the APAR folder, please see https://github.com/nrb171/EXMIRAS for the latest version of this and other scripts.

% use the job scheduler (/glade/u/home/nbarron/workshop/EXMIRAS/scripts/helpers/submitEXMIRASJobs.sh) to run the ideal simulation on Derecho. This should be a fairly good example of how to run this on other HPC systems.
T0 = str2num(getenv('T0'));
humidity = str2num(getenv('HUMIDITY'));
for dBZStart = 30:5:50
    for zdr = 0.4:0.4:2
        parFunction(T0, humidity, dBZStart, zdr, "S", -6.5)
    end
end

%%% function to run the ideal simulation %%%
function parFunction(T0m, humidity, dBZStart, zdr, bandName, lapseRate)
    % Run the ideal simulation                  
    %%! set up initial droplet size distribution
    fprintf('T0m = %1.2f, humidity = %1.2f, dBZStart = %d, zdr = %1.2f\n', T0m, humidity, dBZStart, zdr)
    ex = exmiras;
    %% set up initial droplet size distribution
    ex.xgrid = [50];
    ex.ygrid = [50];
    ex.zgrid = 50:50:2050;

    sx = numel(ex.xgrid);
    sy = numel(ex.ygrid);
    sz = numel(ex.zgrid);
    dBZi = dBZStart + rand(sx, sy, sz)*0;
    Zdri = zdr + rand(sx, sy, sz)*0;
    dBZi(:,:,1:size(dBZi,3)-1) = -inf;
    Zdri(:,:,1:size(dBZi,3)-1) = 0;

    ex = ex.initFromLambdaName("S");
    ex = ex.initFromReflectivity(dBZi, Zdri);

    ex.NP = ex.N;

    %%! set up initial state variables
    temp = @(z) lapseRate*z/1000 + T0m;
    pres = @(z) 1000*exp(-z/1000/8.4);
    es = @(T) 2.53e9*exp(-5420./(T)); %hPa
    [ym, xm, zm] = meshgrid(ex.ygrid, ex.xgrid, ex.zgrid);
    ex.T = temp(zm);
    ex.p = pres(zm);
    ex.pv = es(ex.T)*humidity;

    ex.u = zeros(size(ex.N, [1:3]));
    ex.v = zeros(size(ex.N, [1:3]));
    ex.w = zeros(size(ex.N, [1:3]));

    ex.dt = 0.5;
    ex.st = 5*3600./ex.dt; 

    %% set storm time
    numSteps = ex.st+0;
    variablesToSave = {'T', 'p', 'pv', 'qv', 'Zhh', 'Zvv', 'RR', 'rhohv', 'kdp', 'Zdr', 'theta'};
    ExmirasRun = struct();
    ExmirasRun.initVariables.p = pres(zm);
    ExmirasRun.initVariables.T = temp(zm);
    ExmirasRun.initVariables.pv = es(ex.T)*humidity;
    ExmirasRun.ID = sprintf('%s_ideal_%d_%1.2f_%2.0d_%1.2f', bandName, T0m, humidity,dBZStart, zdr);

    %% check to see if the file already exists, if so, break the loop
    filesInSaveDir = dir('/glade/u/home/nbarron/work/exmirasData/');
    if any(contains({filesInSaveDir.name}, ExmirasRun.ID))
        fprintf('File %s already exists, skipping...\n', ExmirasRun.ID)
        return
    end

    %% set up the struct to save the run
    for i = 1:numel(variablesToSave)
        ExmirasRun.(variablesToSave{i}) = NaN(numSteps, sz);
    end

    %% integrate the model forward
    for i = 1:numSteps
        if mod(i/numSteps, 0.05) == 0
            fprintf('%.2f%%\n', i/numSteps*100)
        end
        ex = ex.integrate;

        % save variables
        for j = 1:numel(variablesToSave)
            ExmirasRun.(variablesToSave{j})(i,:) = squeeze(ex.(variablesToSave{j})(1,1,:));
        end
        ExmirasRun.Zdr(i,:) = ex.Zhh(1,1,:)-ex.Zvv(1,1,:);

        % stop integrating if the water vapor saturation is abov 98%
        if all(ex.qv(1,1,:)>=0.98)
            % remove all unused time steps
            for j = 1:numel(variablesToSave)
                ExmirasRun.(variablesToSave{j})(i+1:end,:) = [];
            end
            ExmirasRun.Zdr(i+1:end,:) = [];

            % end integration
            fprintf('Water vapor saturation reached at time step %d\n', i)
            break
        end 
    end

    % set up grids/save run
    ExmirasRun.Grids.tgrid = ((1:i))*ex.dt;
    ExmirasRun.Grids.zgrid = ex.zgrid;
    save(['/glade/u/home/nbarron/work/exmirasData/',ExmirasRun.ID,'.mat'], 'ExmirasRun')

    % Plot the height-time profiles for the chosen variables
    [TT, ZT] = meshgrid(ExmirasRun.Grids.tgrid, ExmirasRun.Grids.zgrid);

    variablesToPlot = {'T', 'qv', 'Zhh', 'Zdr', 'RR', 'rhohv', 'kdp'};
    titles = {'temperature perturbation', 'water vapor saturation', 'Z_{HH}', 'Z_{DR}', 'Rain Rate', 'Correlation Coefficient', 'Specific Differential Phase'};
    cLabels = {'K', '', 'dBZ', 'dB', 'mm/hr', ' ', 'deg/km'};
    clims = {[270,310], [0.1,1], [-2*10/6,50], [0,8], [0, 8], [0.99, 1], [0.5e-2, 0.5e-1]};
    colormaps = {colortables('newblue'), flipud(winter(64)), colortables('RADAR32'), flipud(hsv(64)), flipud(spring(64)), flipud(turbo(64)), flipud(autumn(64))+[0,0,0.8]};
    for i = 1:numel(variablesToPlot)
        fig = figure("Units", "inches", "Position", [0,0,5.75,1.5]);

        pcolor(TT, ZT, real(ExmirasRun.(variablesToPlot{i}))');
        shading interp
        xlabel('t [s]')
        ylabel('z [m]')
        title(titles{i})
        cb = colorbar;
        cb.Label.String = cLabels{i};
        colormap(colormaps{i})
        clim(clims{i})
        shading flat
        titles{i} = strrep(strrep(titles{i}, "{", ''), "}", '');
        print2(fig, strrep(sprintf('/glade/u/home/nbarron/work/exmirasData/Figures/%s-%s.png', ExmirasRun.ID, titles{i}), ' ', '-'))
    end
end % end of parFunction

