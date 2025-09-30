

folders = dir('/h/eol/nbarron/work/pecan/spol/rhi/2015*');
for di = 4:length(folders)
    d = "/h/eol/nbarron/work/pecan/spol/rhi/"+folders(di).name+"/";
    files = dir(d+"*cfrad*");

    zGrid = 25:50:2025;
    xGrid = 3000:250:100000;

    radarLat = 37;
    radarLon = -98.5;


    % fig1=figure('Units', 'inches', 'Position', [0,0,3,3])
    % ax1 = axes(fig1);
    % hold on
    % muFun = @(gamma) -0.016*gamma.^2 + 1.213*gamma - 1.957;
    % plot(0:20, muFun(0:20), 'k--')

    % muFun2 = @(gamma)-0.0201.*gamma.^2 + 0.902.*gamma -  1.718;
    % plot(0:20, muFun2(0:20), 'r--')

    % fig2 = figure;
    % ax2 = axes(fig2);

    dz=100;
    [zrmesh, rzmesh] = meshgrid(zGrid, xGrid);

    Dmean = 1.1;
    vmean = -0.1021 + 4.932*Dmean - 0.955*Dmean.^2 + 0.07934*Dmean.^3 - 0.002362*Dmean.^4;


    %% build map of sounding times
    fid = fopen("/h/eol/nbarron/work/pecan/upper-air/composite/PECAN_5MB_"+folders(di).name+".cls");
    sID = 1;
    str = fgetl(fid);
    sounding = {};
    
    while ~isempty(str)
        str = fgetl(fid);
        %% look for start of a sounding

        if contains(str, "UTC Release Time (y,m,d,h,m,s):")
            soundingTimes(sID) = datetime(sscanf(str, "UTC Release Time (y,m,d,h,m,s): %f, %f, %f, %f:%f:%f")');
        end
        if contains(str, '------')
            str = fgetl(fid);
            sounding{sID} = NaN(1,20);
            %% loop until the next "Data" line
            while ~contains(str, "Data")
                sounding{sID}(end+1,:) = sscanf(str, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", [20,1])';
                str = fgetl(fid);

                if all(isnan(sounding{sID}(1,:)))
                    sounding{sID}(1,:) = [];
                end
                if feof(fid)
                    break
                end
            end
            sID = sID + 1;
        end
        if feof(fid)
            break
        end
    end

    soundingLats = cellfun(@(x) x(1,12), sounding);
    soundingLons = cellfun(@(x) x(1,11), sounding);

    %% find the distance from soundings to radar
    dists = distance(radarLat, radarLon, soundingLats, soundingLons);

    inds = find(dists < 1)
    try
        Soundings.RHsoundings = reshape(cell2mat(cellfun(@(x) interp1(x(:,15)-x(1,15),x(:,5), 25:50:2025), sounding(inds), 'UniformOutput', false)), [], 41)'; % relative humidity
    catch
        Soundings.RHsoundings = reshape(cell2mat(cellfun(@(x) interp1(x(:,15)-x(1,15)+rand(size(x(:,15)))*0.001,x(:,5), 25:50:2025), sounding(inds), 'UniformOutput', false)), [], 41)'; % relative humidity
        Soundings.RHsoundings(Soundings.RHsoundings >120) = NaN;
        Soundings.RHsoundings = fillmissing(Soundings.RHsoundings, 'pchip', 1);
    end
    Soundings.timesInDatetime = soundingTimes(inds);
    Soundings.dists = dists(inds);
    Soundings.lats = soundingLats(inds);
    Soundings.lons = soundingLons(inds);

    %% set up bin edges and storage for mu/gamma histogram
    gammaEdges = 0:0.25:20;
    muEdges = -2:0.25:15;
    N = zeros(numel(muEdges)-1, numel(gammaEdges)-1);

    gammaRange = linspace(0, 20, 100);

    parfor ii = 1:numel(files)
        
            %% load necessary rhi data
            nci = ncinfo(fullfile(files(ii).folder,files(ii).name));
            time = ncread(fullfile(files(ii).folder,files(ii).name), 'time');
            rayStartIndex = ncread(fullfile(files(ii).folder,files(ii).name), 'ray_start_index');
            rayStartRange = ncread(fullfile(files(ii).folder,files(ii).name), 'ray_start_range');
            rayNGates = ncread(fullfile(files(ii).folder,files(ii).name), 'ray_n_gates');
            range = ncread(fullfile(files(ii).folder,files(ii).name), 'range');
            azimuth = ncread(fullfile(files(ii).folder,files(ii).name), 'azimuth');
            elevation = ncread(fullfile(files(ii).folder,files(ii).name), 'elevation');
            Zhht = ncread(fullfile(files(ii).folder,files(ii).name), 'DBZ_F');
            
            %% remove the file if there's no precipitation
            if max(Zhht) < 30
                % fprintf('File %s has no precipitation, deleting...\n', files(ii).name)
                % delete(fullfile(files(ii).folder,files(ii).name))
                continue
            end

            %% convert to 2D: azimuth vs. gate space
            Zdrt = ncread(fullfile(files(ii).folder,files(ii).name), 'ZDR_F');
            Zhh = NaN(max(rayNGates), numel(rayNGates));
            Zdr = NaN(max(rayNGates), numel(rayNGates));
            for jj = 1:numel(rayNGates)
                Zhh(1:rayNGates(jj),jj) = Zhht(rayStartIndex(jj)+1:rayStartIndex(jj)+rayNGates(jj));
                Zdr(1:rayNGates(jj),jj) = Zdrt(rayStartIndex(jj)+1:rayStartIndex(jj)+rayNGates(jj));
            end



            %% get datetime
            dt = datetime(files(ii).name(7:21), 'InputFormat', 'yyyyMMdd_HHmmss');

            %% interpolate RH profile to time of RHI with exponential weights
            w = exp(-abs(Soundings.dists.*(hours(Soundings.timesInDatetime - dt))/3));
            RHProfile = sum(Soundings.RHsoundings.*w, 2)./sum(w)/100;

            %% set up dsd estimator
            de = dsdAssimilation('S');
            % de.bandName = 'S';
            de.RHProfileObs = RHProfile;    

            %% loop over all azimuths
            noPrecipInAzimuthCounter = 0;
            for az = unique(round(azimuth/5)*5)'
                fprintf('Processing file %s, azimuth %1.0f...\n', files(ii).name, az)
                mask = abs(azimuth-az)<5; % only look at azimuths near 180 deg
                if nnz(mask) < 10
                    noPrecipInAzimuthCounter = noPrecipInAzimuthCounter + 1;
                    continue
                end
                [elevmesh,rangemesh] = meshgrid(elevation(mask), range);
                xmesh = rangemesh.*cosd(elevmesh);
                ymesh = rangemesh.*sind(elevmesh);
                Zhhm = Zhh(:,mask);
                Zdrm = Zdr(:,mask);

                if max(Zhhm(:)) < 30
                    noPrecipInAzimuthCounter = noPrecipInAzimuthCounter + 1;
                    continue

                end

                Zhhsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(Zhhm(:)), 'linear');
                Zdrsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(Zdrm(:)), 'linear');

                ZhhProfile = Zhhsi(rzmesh, zrmesh);
                ZdrProfile = Zdrsi(rzmesh, zrmesh);


                
                %% perform the DSD estimation only when the Zhh is above 30 dBZ
                inds = find(ZhhProfile(:,end)>30 & max(ZhhProfile(:,:), [], 2)<60);
                inds = inds(1:3:end); % only take every third point to speed up the calculation
                if isempty(inds)
                    noPrecipInAzimuthCounter = noPrecipInAzimuthCounter + 1;
                    
                    continue
                end

                
                muTop = [];
                lambdaTop = [];
                N0Top = [];
                error = [];
                rhmean = [];
                ZhhTop = [];
                ZdrTop = [];

                muBot = [];
                lambdaBot = [];
                N0Bot = [];
                ZhhBot = [];
                ZdrBot = [];

                gateRange = [];
                ra = radar('S');
                for jj = 1:numel(inds)
                    indToTest = max(inds(jj)-1,1): min(inds(jj)+1, size(ZhhProfile,1));
                    
                    [N0TopTemp, muTopTemp, lambdaTopTemp,fv]=de.profileOptimizer(...
                        ZhhProfile(indToTest,:)', ... 
                        ZdrProfile(indToTest,:)' ...
                    );

                    [dN,~,~] = de.estimateSingleRadarProfile(N0TopTemp, muTopTemp, lambdaTopTemp);

                    
                    
                    % keyboard
                    
                    [N0BotTemp, muBotTemp, lambdaBotTemp] = de.getN0MuLambdaFromN(dN(1,:));
                    % plot(de.D, dN(1,:))
                    % hold on
                    % plot(de.D, de.getNFromN0MuLambda(N0Bot, muBot, lambdaBot), 'r--')
                    % legend('Top', 'Bottom')
                    % hold off
                    % print2
                    % ra.rngToggle = true; % toggle for random number generator
                    % ra.rngToggle = false; % toggle for random number generator
                    % [~, ~,muTemp, gammaTemp] = ra.initFromDm(ZhhOpt, ZdrOpt, DmOpt);
                    muTop(jj) = muTopTemp;
                    lambdaTop(jj) = lambdaTopTemp;
                    N0Top(jj) = N0TopTemp;
                    ZhhTop(jj) = ra.calcZhh(dN(41,:));
                    ZdrTop(jj) = ra.calcZdr(dN(41,:));
                    error(jj) = fv;

                    muBot(jj) = muBotTemp;
                    lambdaBot(jj) = lambdaBotTemp;
                    N0Bot(jj) = N0BotTemp;
                    ZhhBot(jj) = ra.calcZhh(dN(1,:));
                    ZdrBot(jj) = ra.calcZdr(dN(1,:));
                    
                    rhmean(jj) = mean(RHProfile);
                    gateRange(jj) = rzmesh(inds(jj), end);


                end

                %% create mu/gamma database
                s = struct('ZhhTop', ZhhTop, 'ZdrTop', ZdrTop, 'muTop', muTop, 'lambdaTop', lambdaTop, 'N0Top', N0Top, ...
                    'ZhhBot', ZhhBot, 'ZdrBot', ZdrBot, 'muBot', muBot, 'lambdaBot', lambdaBot, 'N0Bot', N0Bot, ...
                    'error', error, 'rhmean', rhmean, 'gateRange', gateRange);

                save(sprintf('/h/eol/nbarron/work/dsd-estimator/spol-N0MuLambda-%s-%1.0f.mat', files(ii).name(7:21), az), '-fromstruct', s)
                
                % scatter(ax1,gamma, mu)
                % title(ax1, files(ii).name)
                % print2(fig1, ['~/figures/publications/exmiras/muvlambda-scatter-', folders(di).name, '.png'])
                
                % N = N+histcounts2(gamma,mu,gammaEdges, muEdges)';
                
            end
            if noPrecipInAzimuthCounter == numel(unique(round(azimuth/5)*5)')
                fprintf('File %s has no precipitation, deleting...\n', files(ii).name)
                delete(fullfile(files(ii).folder,files(ii).name))
            end

    fprintf('Finished file %s\n', files(ii).name)
    end
end
quit


%% load all of the mu/gamma data and make a plot
%%! need to filter by rhmean
    files = dir('/h/eol/nbarron/work/dsd-estimator/spol-N0MuLambda*.mat');
    muTop = [];
    lambdaTop = [];
    N0Top = [];
    ZhhTop = [];
    ZdrTop = [];
    muBot = [];
    lambdaBot = [];
    N0Bot = [];
    ZhhBot = [];
    ZdrBot = [];
    error = [];
    rhmean = [];
    gateRange = [];

    gammaEdges = 0:0.25:20;
    gammaRange = linspace(0, 20, 100);
    muEdges = -2:0.25:15;
    N = zeros(numel(muEdges)-1, numel(gammaEdges)-1);

    for ii = 1:numel(files)
        try
            f = load(fullfile(files(ii).folder, files(ii).name));

            if nnz(fieldnames(f) == ["ZhhTop", "ZdrTop", "muTop", "lambdaTop", "N0Top", ...
                    "ZhhBot", "ZdrBot", "muBot", "lambdaBot", "N0Bot", ...
                    "error", "rhmean", "gateRange"]) ~= 13
                fprintf('File %s is missing some fields, skipping...\n', files(ii).name)
                continue
            end
            
            muTop = [muTop, f.muTop];
            lambdaTop = [lambdaTop, f.lambdaTop];
            N0Top = [N0Top, f.N0Top];
            ZhhTop = [ZhhTop, f.ZhhTop];
            ZdrTop = [ZdrTop, f.ZdrTop];
            muBot = [muBot, f.muBot];
            lambdaBot = [lambdaBot, f.lambdaBot];
            N0Bot = [N0Bot, f.N0Bot];
            ZhhBot = [ZhhBot, f.ZhhBot];
            ZdrBot = [ZdrBot, f.ZdrBot];
            error = [error, f.error];
            rhmean = [rhmean, f.rhmean];
            gateRange = [gateRange, f.gateRange];

        end
    end

    figure('Units', 'inches', 'Position', [0,0,3,3])
    % [gammaMesh, muMesh] = meshgrid((gammaEdges(1:end-1) + gammaEdges(2:end))/2, (muEdges(1:end-1) + muEdges(2:end))/2);
    % gammaRescaled = [];
    % muRescaled = [];
    % for jj = 1:numel(gammaMesh)
    %     N2 = round(log10(1+N(jj)));
    %     if N2 == 0 
    %         continue
    %     end
    %     gammaRescaled(end+1:end+N2) = gammaMesh(jj);
    %     muRescaled(end+1:end+N2) = muMesh(jj);
    % end
    for rh = 0.4:0.2:0.9
        mask = ZhhTop > 30 & ZhhTop < 100 & rhmean > rh & rhmean <rh+0.2;
        N = histcounts2(lambdaTop(mask), muTop(mask), gammaEdges, muEdges)';
        
        figure('Units', 'inches', 'Position', [0,0,3,3])
        hold off
        pcolor((gammaEdges(1:end-1) + gammaEdges(2:end))/2, (muEdges(1:end-1) + muEdges(2:end))/2, log10(N), 'EdgeColor', 'none')
        hold on
        muFun = @(gamma) -0.016*gamma.^2 + 1.213*gamma - 1.957;
        plot(gammaRange, muFun(gammaRange), 'k--')
        gammaRange = linspace(0, 20, 100);
        title("top")
        muFun2 = @(gamma)-0.0201.*gamma.^2 + 0.902.*gamma -  1.718;
        plot(gammaRange, muFun2(gammaRange), 'r--')
        print2
        keyboard
        
        
        NBot = histcounts2(lambdaBot(mask), muBot(mask), gammaEdges, muEdges)';
        figure('Units', 'inches', 'Position', [0,0,3,3])
        hold off
        pcolor((gammaEdges(1:end-1) + gammaEdges(2:end))/2, (muEdges(1:end-1) + muEdges(2:end))/2, log10(NBot), 'EdgeColor', 'none')
        hold on
        plot(gammaRange, muFun(gammaRange), 'k--')
        plot(gammaRange, muFun2(gammaRange), 'r--')
        title("bot")
        print2
        keyboard

        % gammaRescaled(muRescaled < 0) = [];
        % muRescaled(muRescaled < 0) = [];
        % muRescaled(gammaRescaled < 1) = [];
        % gammaRescaled(gammaRescaled < 1) = [];

        % pf = polyfit(gamma(mask), mu(mask), 2)
        % errorfun = @(x) abs((mean(nonzeros(((polyval(x, gammaMesh) - muMesh).^2).*log10(1+N)))).^(1/2));
        % errorPoly = rms((polyval(pf, gammaMesh) - muMesh).*(N), "all")
        % pf=fminsearchbnd(errorfun, [0,0,0], [-0.25, -3, -3], [0.25, 3, 3])
        % % pf=fminsearchbnd(errorfun, [0,0,0])
        % plot(gammaRange, polyval(pf, gammaRange), 'Color',[94,79,162]/255, 'LineWidth', 1.5)
        % xlim([0,20])
        % ylim([-2,14])   
        % colormap(cbrewer('seq', 'YlGnBu', 256))

        % otherFit = fit(gamma(mask)', mu(mask)', 'power2');
        % errorPowerLaw = rms((otherFit(gammaMesh) - muMesh(:)).*(N(:)), "all")
        % plot(gammaRange, otherFit(gammaRange), 'Color',[0.9290, 0.6940, 0.1250], 'LineWidth', 1.5)
        % print2(sprintf('~/figures/publications/exmiras/dsdEstimatorMuvGamma-spol-%0.1f.pdf', rh))
        % xlabel('$\Gamma$')
        % ylabel('$\mu$')
        % print2(gcf, './.temp3.png')
    end


    X = [Zhh; Zdr; gamma; mu; rhmean]';
    varNames = ["Z_{hh} (dBZ)", "Z_{dr} (dB)", "\Lambda", "\mu", "RH_{mean}"];

    mask = error < 0.2 & Zhh > 40 & Zhh < 50 & gateRange < 80000;
    rhclass = discretize(rhmean, 0.4:0.2:1);
    figure('Units', 'inches', 'Position', [0,0,8,8])
    gplotmatrix(X(mask',:),[],rhclass(mask),[],[],8,[],[],varNames)
    print2




    figure('Units', 'inches', 'Position', [0,0,0.25,3])
    colormap(cbrewer('seq', 'YlGnBu', 256))
    ax = gca;
    ax.Visible = 'off';
    cb = colorbar;
    clim([1+min(N(:)), 1+max(N(:))])
    set(gca,'ColorScale','log')
    cb.Ticks = [10,50,100,250,500,850];
    print2('~/figures/publications/exmiras/dsdEstimatorMuvGammaCB.pdf')

