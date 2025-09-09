
%% generate simulations to linearize the DSD evolution
% this script generates a series of idealized EXMIRAS simulations to
% characterize the evolution of the DSD as a function of initial Dm, Zdr,
% and RH. The reflectivity is solved in the forward model in `dsdAssimilation.m`

% To generate the simulations, run the following loop. This will take about 4 weeks at this resolution.

% To save the LUTs, run the code beginning at line 34.
zgrid = 25:50:2025;
lapseRate = -6.5; % K/km

for Dm = 0.1:0.1:1.8
    for zdr = 0.3:0.3:3
        for RH = 0.1:0.05:0.95
            for bandName = {'S', 'C', 'X'}
                temp = @(z) lapseRate*z/1000 + 290;
                pres = @(z) 1000*exp(-z/1000/8.4);

                bandName = bandName{1};
                % T0m = 290; % K
                saveDir = '~/work/exmirasData/LUTs-';
                T = temp(zgrid);
                p = pres(zgrid);
                runEXMIRAS(40, zdr, Dm, bandName, temp(zgrid), pres(zgrid), ones(size(zgrid))*RH, saveDir, 'waterVaporInteraction', false)
            end
           
        end
    end
end

return

%% generate the LUTs
ZhhStart = 30;

% load all of the files
files = dir('~/work/exmirasData/LUTs-*');
files = {files.name};

% loop over the minutes of simulation.
for midx = [4,5,7,10,15]
    tic

    %% assign grids for the interpolant
    zdrgrid = round(0.3:0.3:3,1);
    Dmgrid = round(0.1:0.1:1.8,1);
    RHgrid = round(0.1:0.05:0.95,2);
    Ntop3 = NaN(numel(zdrgrid), numel(RHgrid), numel(Dmgrid), 250);
    Nbot3 = NaN(numel(zdrgrid), numel(RHgrid), numel(Dmgrid), 250);

    deta =struct("S", NaN(numel(zdrgrid), numel(RHgrid), numel(Dmgrid), 250), ...
                "C", NaN(numel(zdrgrid), numel(RHgrid), numel(Dmgrid), 250), ...
                "X", NaN(numel(zdrgrid), numel(RHgrid), numel(Dmgrid), 250));

    % loop over S-, C-, and X-band
    for bandName = {'S', 'C', 'X'}
        filesBand = files(contains(files, ['LUTs-', bandName{1}]));
        for fname=filesBand
            fileInfo.Dm = round(str2num(fname{1}(14:17)),1);
            fileInfo.zdr = round(str2num(fname{1}(19:22)),1);
            fileInfo.RH = round(str2num(fname{1}(24:27)),1);

            if ~all([...
                any(zdrgrid == fileInfo.zdr), ...
                any(RHgrid == fileInfo.RH), ...
                any(Dmgrid == fileInfo.Dm)])

                continue
            end
            f = load(['~/work/exmirasData/', fname{1}]);

            % calculate deta as (Nbot-Ntop)/dz/Ntop
            Ntop = squeeze(f.ExmirasRun.N(midx,end,:));
            Nbot = squeeze(f.ExmirasRun.N(midx,1,:));
            Nbot(1) = Nbot(2);
            deta.(bandName{1})(...
                zdrgrid == fileInfo.zdr, ...
                RHgrid == fileInfo.RH, ...
                Dmgrid == fileInfo.Dm, ...
            :) = (Nbot-Ntop)./2000/max(Ntop);

            
            Nbot3(...
                zdrgrid == fileInfo.zdr, ...
                RHgrid == fileInfo.RH, ...
                Dmgrid == fileInfo.Dm, ...
                :) = Nbot;
        end

        
    end
    fprintf('Finished midx = %d\n', midx)
    save(sprintf('../data/LUTs/detaLUTs-hires-%1.0f.mat', midx), 'deta','zdrgrid', 'RHgrid', 'Dmgrid')
    toc
end

%%% function to run the ideal simulation %%%
function runEXMIRAS(dBZStart, zdr, Dm, bandName, T, p, RH, saveDir, varargin)
    tic
    % Run the ideal simulation                  
    %%! set up initial droplet size distribution
    % fprintf('T0m = %1.2f, humidity = %1.2f, dBZStart = %d, zdr = %1.2f\n', T0m, humidity, dBZStart, zdr)
    ex = exmiras;
    ex.rngToggle = false;
    %% set up initial droplet size distribution
    ex.xgrid = [50];
    ex.ygrid = [50];
    ex.zgrid = 25:50:2025;

    sx = numel(ex.xgrid);
    sy = numel(ex.ygrid);
    sz = numel(ex.zgrid);
    dBZi = dBZStart + rand(sx, sy, sz)*0;
    Zdri = zdr + rand(sx, sy, sz)*0;
    dBZi(:,:,1:size(dBZi,3)-1) = -inf;
    Zdri(:,:,1:size(dBZi,3)-1) = 0;

    ex = ex.initFromLambdaName(bandName);
    % ex = ex.initFromReflectivity(dBZi, Zdri);
    ex = ex.initFromDm(dBZStart, zdr, Dm);
    ex.Zdr(end)
    

    ex.NP = ex.N;

    %%! set up initial state variables
    % temp = @(z) lapseRate*z/1000 + T0m;
    % pres = @(z) 1000*exp(-z/1000/8.4);
    es = @(T) 2.53e9*exp(-5420./(T)); %hPa
    % [ym, xm, zm] = meshgrid(ex.ygrid, ex.xgrid, ex.zgrid);
    ex.T = T;
    ex.p = p;
    ex.pv = es(ex.T).*RH;

    ex.T = reshape(ex.T, [1,1,numel(ex.zgrid)]);
    ex.p = reshape(ex.p, [1,1,numel(ex.zgrid)]);
    ex.pv = reshape(ex.pv, [1,1,numel(ex.zgrid)]);

    ex.u = zeros(size(ex.N, [1:3]));
    ex.v = zeros(size(ex.N, [1:3]));
    ex.w = zeros(size(ex.N, [1:3]));

    ex.dt = 0.5;
    % ex.st = 3600./ex.dt; 
    ex.st = 1200./ex.dt; 
    
    % keyboard


    %% set storm time
    numSteps = ex.st+0;
    variablesToSave = {'T', 'p', 'pv', 'qv', 'Zhh', 'Zvv', 'RR', 'rhohv', 'kdp', 'Zdr', 'theta'};
    ExmirasRun = struct();
    ExmirasRun.initVariables.p = p;
    ExmirasRun.initVariables.T = T;
    ExmirasRun.initVariables.pv = es(ex.T).*RH;
    ExmirasRun.ID = sprintf('%s_%s_%s_%s_%s', bandName, "ideal", num2str(Dm, '%1.2f'), num2str(zdr, '%1.2f'), num2str(mean(RH), '%1.2f'));

    fprintf(ExmirasRun.ID)

    %% check to see if the file already exists, if so, break the loop
    % '~/work/exmirasData/LUTs-'
    sdsplit = strsplit(saveDir, '/');
    filesInSaveDir = dir([strjoin(sdsplit(1:end-1),'/'), '/']);
    % keyboard
    if any(contains({filesInSaveDir.name}, ExmirasRun.ID))
        fprintf('File %s already exists, skipping...\n', ExmirasRun.ID)
        return
    end

    %% set up the struct to save the run
    for i = 1:numel(variablesToSave)
        ExmirasRun.(variablesToSave{i}) = NaN(numSteps, sz);
    end
    ExmirasRun.N = [];
    ExmirasRun.dNevap = [];
    ExmirasRun.dNcoal = [];
    ExmirasRun.dNfallout = [];
    ExmirasRun.gamma = [];
    ExmirasRun.mu = [];
    ExmirasRun.N0 = [];

    if numel(varargin)>1
        for ii = 1:numel(varargin)/2
            ex.((varargin{2*ii-1})) = varargin{2*ii};
        end
    end

    % keyboard
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

        if mod(i, 120) == 0
            fprintf('saving dsd at time step %d\n', i)
            ExmirasRun.N(end+1,:,:) = squeeze(ex.N(1,1,:,:));
            ExmirasRun.dNevap(end+1,:,:) = squeeze(ex.dNevap(1,1,:,:));
            ExmirasRun.dNcoal(end+1,:,:) = squeeze(ex.dNcoal(1,1,:,:));
            ExmirasRun.dNfallout(end+1,:,:) = squeeze(ex.dNfallout(1,1,:,:));
            ExmirasRun.gamma(end+1) = ex.gamma(1,1,end);
            ExmirasRun.mu(end+1) = ex.mu(1,1,end);
            ExmirasRun.N0(end+1) = ex.N00(1,1,end);
            % ExmirasRun.dN = squeeze(ex.dN(1,1,:,:
        end

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
    save([saveDir,ExmirasRun.ID,'.mat'], 'ExmirasRun')
    % load([saveDir,'S_ideal_300_0.30_40_0.50','.mat'], 'ExmirasRun')

    % % Plot the height-time profiles for the chosen variables
    % [TT, ZT] = meshgrid(ExmirasRun.Grids.tgrid, ExmirasRun.Grids.zgrid);

    % [~,dTdt] = gradient(ExmirasRun.T, 0.5); % convert to K/hr
    % ExmirasRun.dTdt = dTdt*3600;
    % ExmirasRun.RH = ExmirasRun.qv*100;
    % variablesToPlot = {'dTdt', 'RH', 'Zhh', 'Zdr', 'RR', 'rhohv', 'kdp'};
    % titles = {'CoolingRate', 'water vapor saturation', 'Z_{HH}', 'Z_{DR}', 'Rain Rate', 'Correlation Coefficient', 'Specific Differential Phase'};

    % plotTimeProfile(ExmirasRun.Grids.tgrid, ExmirasRun.Grids.zgrid, ['~/figures/publications/exmiras/',ExmirasRun.ID], ...
    %     cellfun(@(x) ExmirasRun.(x), variablesToPlot, 'uniform', 0), variablesToPlot)
    % cLabels = {'K', '', 'dBZ', 'dB', 'mm/hr', ' ', 'deg/km'};
    % clims = {[270,310], [0.1,1], [-2*10/6,50], [0,8], [0, 8], [0.99, 1], [0.5e-2, 0.5e-1]};
    % colormaps = {colortables('newblue'), flipud(winter(64)), colortables('RADAR32'), flipud(hsv(64)), flipud(spring(64)), flipud(turbo(64)), flipud(autumn(64))+[0,0,0.8]};
    % for i = 1:numel(variablesToPlot)
        % fig = figure("Units", "inches", "Position", [0,0,6.05,1.5]);

        % nuimagesc(TT, ZT, real(ExmirasRun.(variablesToPlot{i}))');
        % isc=imagesc(gca, ZT, TT, real(variablesToPlot{i}'));
        % isc.AlphaData = ones(size(isc.CData));
        % isc.AlphaData(isnan(real(variablesToPlot{iPlot}'))) = 0;
        % shading interp
        % xlabel('t [s]')
        % ylabel('z [m]')
        % title(titles{i})
        % % cb = colorbar;
        % setMeteoColormap(gca, variablesToPlot{i})
        % cb.Label.String = cLabels{i};
        % colormap(colormaps{i})
        % clim(clims{i})
        % shading flat
        % titles{i} = strrep(strrep(titles{i}, "{", ''), "}", '');
        % print2(fig, strrep(sprintf('~/figures/publications/exmiras/%s-%s.pdf', ExmirasRun.ID, variablesToPlot{i}), ' ', '-'))
    % end
    toc
end % end of runEXMIRAS

