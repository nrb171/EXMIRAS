
%% generate simulations to linearize the DSD evolution
% this script generates a series of idealized EXMIRAS simulations to
% characterize the evolution of the DSD as a function of initial mu, lambda,
% and Zhh. The reflectivity profiles are solved in the forward model from `dsdAssimilation.m`
% To generate the simulations, run the following loop. 
% To generate the LUTs from the results of the simulations, run the code from L47--L93.

zgrid = 25:100:2025;
lapseRate = -6.5; % K/km


% mu vs. Lambda bounds.
polytop = [2.3277, -0.9189];
polybot = [0.2796, -2.0];

numSteps = 20;

Zhh = 30;
da = dsdAssimilation('S');

for Zhh = 30%[10, 20, 30:5:60]
    for lambda = linspace(0,20, numSteps)
        muTop = min(polyval(polytop, lambda), 15);
        muBot = polyval(polybot, lambda);
        muRange = linspace(muBot, muTop, numSteps);
        for mu = muRange
            N0 = da.getN0FromZhhMuLambda(Zhh, mu, lambda);
            for RH = 0.2:0.2:1

                % set temperature and pressure profiles
                temp = @(z) lapseRate*z/1000 + 295;
                pres = @(z) 1000*exp(-z/1000/8.4);
                
                % T0m = 290; % K
                saveDir = '~/work/exmirasData/LUTs-hires-N0MuLambda-';
                T = temp(zgrid);
                p = pres(zgrid);
                runEXMIRAS(N0, mu, lambda,  temp(zgrid), pres(zgrid), ones(size(zgrid))*RH, saveDir, 'waterVaporInteraction', false)

            
            end
        end
    end
end

return

%% generate the LUTs
ZhhStart = 30;

% load all of the files
files = dir('~/work/exmirasData/LUTs-hires-N0MuLambda-*');
files = {files.name};

% loop over the minutes of simulation.
for midx = [4,5,7,10,15]
    N0 = [];
    mu = [];
    lambda = [];
    RH = [];
    deta = [];
    tic

    %% assign grids for the interpolant
    zdrgrid = round(0.3:0.3:3,1);
    Dmgrid = round(0.1:0.1:1.8,1);
    RHgrid = round(0.1:0.05:0.95,2);
    Ntop3 = NaN(numel(zdrgrid), numel(RHgrid), numel(Dmgrid), 250);
    Nbot3 = NaN(numel(zdrgrid), numel(RHgrid), numel(Dmgrid), 250);

    for fname=files
        
        f = load(['~/work/exmirasData/', fname{1}]);

        N0(end+1) = f.ExmirasRun.initVariables.N0;
        mu(end+1) = f.ExmirasRun.initVariables.mu;
        lambda(end+1) = f.ExmirasRun.initVariables.lambda;
        RH(end+1) = mean(f.ExmirasRun.initVariables.RH);

        % calculate deta as (Nbot-Ntop)/dz/Ntop
        Ntop = squeeze(f.ExmirasRun.N(midx,end,:));
        Nbot = squeeze(f.ExmirasRun.N(midx,1,:));
        Nbot(1) = Nbot(2);
        deta(end+1,:) = (Nbot-Ntop)./2000/max(Ntop);

    end

    fprintf('Finished midx = %d\n', midx)
    detaInterpolant = scatteredInterpolant(mu', lambda', RH', deta, 'linear', 'linear');
    save(sprintf('../data/LUTs/detaLUTs-N0MuLambda-hires-%1.0f.mat', midx), 'deta','N0', 'mu', 'lambda', 'RH', "detaInterpolant")
    toc
end

%%% function to run the ideal simulation %%%
function runEXMIRAS(N0, mu, lambda, T, p, RH, saveDir, varargin)
    tic
    % Run the ideal simulation                  
    %%! set up initial droplet size distribution
    % fprintf('T0m = %1.2f, humidity = %1.2f, dBZStart = %d, zdr = %1.2f\n', T0m, humidity, dBZStart, zdr)
    ex = exmiras;
    ex.rngToggle = false;
    %% set up initial droplet size distribution
    ex.xgrid = [50];
    ex.ygrid = [50];
    ex.zgrid = 25:100:2025;

    sx = numel(ex.xgrid);
    sy = numel(ex.ygrid);
    sz = numel(ex.zgrid);
    
    ex.N = zeros(sx, sy, sz, ex.nBins);
    ex = ex.initializeRadarSimulator('S');
    ex.N(1,1,end,:) = ex.da.getNFromN0MuLambda(N0, mu, lambda);
    % ex.Zdr(end)
    

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
    ex.st = 900./ex.dt; 
    
    % keyboard


    %% set storm time
    numSteps = ex.st+0;
    variablesToSave = {'T', 'p', 'pv', 'qv', 'RR', 'rhohv', 'kdp', 'theta'};
    ExmirasRun = struct();
    ExmirasRun.initVariables.p = p;
    ExmirasRun.initVariables.T = T;
    ExmirasRun.initVariables.pv = es(ex.T).*RH;
    ExmirasRun.initVariables.N0 = N0;
    ExmirasRun.initVariables.mu = mu;
    ExmirasRun.initVariables.lambda = lambda;
    ExmirasRun.initVariables.RH = RH;
    ExmirasRun.ID = sprintf('%1.3f_%1.3f_%1.3f_%1.3f', N0, mu, lambda, mean(RH));

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
        % ExmirasRun.Zdr(i,:) = ex.Zhh(1,1,:)-ex.Zvv(1,1,:);

        if mod(i, 120) == 0
            fprintf('saving dsd at time step %d\n', i)
            ExmirasRun.N(end+1,:,:) = squeeze(ex.N(1,1,:,:));
            ExmirasRun.dNevap(end+1,:,:) = squeeze(ex.dNevap(1,1,:,:));
            ExmirasRun.dNcoal(end+1,:,:) = squeeze(ex.dNcoal(1,1,:,:));
            ExmirasRun.dNfallout(end+1,:,:) = squeeze(ex.dNfallout(1,1,:,:));
            % keyboard
            % ExmirasRun.gamma(end+1) = ex.gamma(1,1,end);
            % ExmirasRun.mu(end+1) = ex.mu(1,1,end);
            % ExmirasRun.N0(end+1) = ex.N00(1,1,end);
            % ExmirasRun.dN = squeeze(ex.dN(1,1,:,:
        end

        % stop integrating if the water vapor saturation is abov 98%
        % if all(ex.qv(1,1,:)>=0.98)
        %     % remove all unused time steps
        %     for j = 1:numel(variablesToSave)
        %         ExmirasRun.(variablesToSave{j})(i+1:end,:) = [];
        %     end
        %     % ExmirasRun.Zdr(i+1:end,:) = [];

        %     % end integration
        %     fprintf('Water vapor saturation reached at time step %d\n', i)
        %     break
        % end 
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

