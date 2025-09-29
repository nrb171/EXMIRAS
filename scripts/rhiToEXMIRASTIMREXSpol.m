%% this script prepares the rhi data for EXMIRAS run, and runs the model
% specifically for the TIMREX-SPOL on 20080614, 

%% Prepare rhi data for EXMIRAS run
    d = "/h/eol/nbarron/work/timrex/spol/20080614/";
    files = dir(d+"*RHI*");

    zGrid = 25:50:2025;
    xGrid = 0:150:150000;

    dz=100;
    [zrmesh, rzmesh] = meshgrid(zGrid, xGrid);

    Dmean = 1.1;
    vmean = -0.1021 + 4.932*Dmean - 0.955*Dmean.^2 + 0.07934*Dmean.^3 - 0.002362*Dmean.^4;


%% get sounding data
    soundingName = "/h/eol/nbarron/work/timrex/upper-air/spol.200806140600.txt";
    ua = readmatrix(soundingName, 'NumHeaderLines', 1, 'FileType', 'text');
    p = ua(:,1); % hPa
    T = ua(:,3); % C
    dwp = ua(:,4); % C


    RH = 100 - 5*(T-dwp);

    TProfile = interp1(ua(:,2)-ua(1,2), T, zGrid, 'linear', 'extrap') + 273.15; % K
    
    RHProfile = interp1(ua(:,2)-ua(1,2), RH, zGrid, 'linear', 'extrap')/100; % relative humidity
    presProfile = interp1(ua(:,2)-ua(1,2), p, zGrid, 'linear', 'extrap'); % pressure hPa
    dwpProfile = interp1(ua(:,2)-ua(1,2), dwp, zGrid, 'linear', 'extrap') + 273.15; % K

    % plot the sounding
    figure("Units","inches","Position",[0,0,3,3])
    skewt(ua(:,1), ua(:,3), ua(:,4))
    print2
    print2(gcf, '/h/eol/nbarron/figures/publications/exmiras/TIMREX_Sounding.pdf')
    
    
    % TProfile = [fliplr(9.8*zGrid(zGrid<ua(1,2)))/1000 + TProfile(1), TProfile]; % extend the profile to the ground
    % RHProfile = [ones(1,numel(zGrid(zGrid<ua(1,2))))*RHProfile(1), RHProfile]; % extend the profile to the ground

    
%% load necessary rhi data
    ii=3
    % ii=9

    nci = ncinfo(fullfile(files(ii).folder,files(ii).name));
    time = ncread(fullfile(files(ii).folder,files(ii).name), 'time');
    rayStartIndex = ncread(fullfile(files(ii).folder,files(ii).name), 'ray_start_index');
    rayStartRange = ncread(fullfile(files(ii).folder,files(ii).name), 'ray_start_range');
    rayNGates = ncread(fullfile(files(ii).folder,files(ii).name), 'ray_n_gates');

    VRt = ncread(fullfile(files(ii).folder,files(ii).name), 'VR');
    Kdpt = ncread(fullfile(files(ii).folder,files(ii).name), 'NKDP');
    range = ncread(fullfile(files(ii).folder,files(ii).name), 'range');
    azimuth = ncread(fullfile(files(ii).folder,files(ii).name), 'azimuth');
    unique(azimuth)
    elevation = ncread(fullfile(files(ii).folder,files(ii).name), 'elevation');
    Zhht = ncread(fullfile(files(ii).folder,files(ii).name), 'DBZ');

    Zdrt = ncread(fullfile(files(ii).folder,files(ii).name), 'ZDR');
    noise = ncread(fullfile(files(ii).folder,files(ii).name), 'unambiguous_range')
    Zhh = NaN(max(rayNGates), numel(rayNGates));
    Zdr = NaN(max(rayNGates), numel(rayNGates));
    VR = NaN(max(rayNGates), numel(rayNGates));
    Kdp = NaN(max(rayNGates), numel(rayNGates));

    for jj = 1:numel(rayNGates)
        Zhh(1:rayNGates(jj),jj) = Zhht(rayStartIndex(jj)+1:rayStartIndex(jj)+rayNGates(jj));
        Zdr(1:rayNGates(jj),jj) = Zdrt(rayStartIndex(jj)+1:rayStartIndex(jj)+rayNGates(jj));
        VR(1:rayNGates(jj),jj) = VRt(rayStartIndex(jj)+1:rayStartIndex(jj)+rayNGates(jj));
        Kdp(1:rayNGates(jj),jj) = Kdpt(rayStartIndex(jj)+1:rayStartIndex(jj)+rayNGates(jj));
    end

%% get datetime and relative humidity profile
    dt = datetime(files(ii).name(7:21), 'InputFormat', 'yyyyMMdd_HHmmss');

    de = dsdAssimilation('S');
    % de.bandName = 'S';
    de.RHProfileObs = RHProfile;
    de.hires = false;

%% extract the azimuth we need
    az = 36
    % az=

    mask = abs(azimuth-az)<0.5; % only look at azimuths near the target
    [elevmesh,rangemesh] = meshgrid(elevation(mask), range);
    xmesh = rangemesh.*cosd(elevmesh);
    ymesh = rangemesh.*sind(elevmesh);
    Zhhm = Zhh(:,mask);
    Zdrm = Zdr(:,mask);
    VRm = VR(:,mask);
    Kdpm = Kdp(:,mask);

%% plot the RHI observations
    meanVel = mean(VRm(Zhhm>30 & VRm>-26.7));

    figure('units','inches','position',[0,0,6.05,1.5])
    pcolor(xmesh./meanVel, ymesh, Zhhm)
    setMeteoColormap(gca, 'Zhh')

    xticks(0:1000:6000)
    xticklabels(round((0:1000:6000)*meanVel))
    shading flat
    ylim([0,2500])
    xlim([-0, 4500])
    print2(gcf, '/h/eol/nbarron/figures/publications/exmiras/TIMREX_Obs_Zhh.pdf')
    figure('units','inches','position',[0,0,6.05,1.5])
    
    
    pcolor(xmesh./meanVel, ymesh, Zdrm)
    setMeteoColormap(gca, 'Zdr')
    shading flat
    ylim([0,2500])

    xlim([-0, 4500])
    xticks(0:1000:6000)
    xticklabels(round((0:1000:6000)*meanVel))
    print2(gcf, '/h/eol/nbarron/figures/publications/exmiras/TIMREX_Obs_Zdr.pdf')

    pcolor(xmesh./meanVel, ymesh, VRm)
    setMeteoColormap(gca, 'RadialVelocity')
    shading flat
    ylim([0,2500])
    xlim([-0, 4500])
    print2(gcf, '/h/eol/nbarron/figures/publications/exmiras/TIMREX_Obs_VR.pdf')

    pcolor(xmesh./meanVel, ymesh, Kdpm)
    setMeteoColormap(gca, 'Kdp')
    shading flat
    ylim([0,2500])
    xlim([-0, 4500])
    xticks(0:1000:6000)
    xticklabels(round((0:1000:6000)*meanVel))
    print2(gcf, '/h/eol/nbarron/figures/publications/exmiras/TIMREX_Obs_Kdp.pdf')

%% load and plot SUR obs.
%{
    folder = "/h/eol/nbarron/work/timrex/spol/20080614/";
    filename = "cfrad.20080614_083002.000_to_20080614_083138.000_SPOLRVP8_v25_SUR.nc";
    

    azimuths = ncread(fullfile(folder,filename), 'azimuth');
    elevations = ncread(fullfile(folder,filename), 'elevation');
    ranges = ncread(fullfile(folder,filename), 'range');
    rayStartIndex = ncread(fullfile(folder,filename), 'ray_start_index');
    rayStartRange = ncread(fullfile(folder,filename), 'ray_start_range');
    rayNGates = ncread(fullfile(folder,filename), 'ray_n_gates');

    ZhhFlat = ncread(fullfile(folder,filename), 'DBZ');

    clear ZhhSur rangeSur azimuthSur elevationSur
    for jj = 1:numel(rayNGates)
        ZhhSur(1:rayNGates(jj),jj) = ZhhFlat(rayStartIndex(jj)+1:rayStartIndex(jj)+rayNGates(jj));
        rangeSur(1:rayNGates(jj),jj) = range(1:rayNGates(jj));
        azimuthSur(1:rayNGates(jj),jj) = azimuths(jj);
        elevationSur(1:rayNGates(jj),jj) = elevations(jj);
    end

    mask = abs(wrapTo360(azimuthSur)-315)<45 & rangeSur<50000;
    [x,y,z] = sph2cart(azimuthSur(mask),elevationSur(mask),rangeSur(mask));
    x = double(x);
    y = double(y);
    z = double(z);

    SI = scatteredInterpolant(x(:), y(:), z(:), reshape(10.^(ZhhSur(mask)/10), [], 1),'linear');

    % directionOfLine = 110
    % xmesh2 = cosd(directionOfLine)*(5000:500:60000) + reshape((-20000:1000:20000).*abs(sind(directionOfLine)),[],1);
    % ymesh2 = sind(directionOfLine)*(5000:500:60000) + reshape((-20000:1000:2000).*abs(cosd(directionOfLine)),[],1);
    % zmesh2 = 0:50:2000;
    % xmesh3 = repmat(xmesh2, 1, 1,numel(zmesh2));
    % ymesh3 = repmat(ymesh2, 1, 1,numel(zmesh2));
    % % zmesh3 = repmat(reshape(zmesh2,1,1,[]
    % [~, ~, zmesh3] = meshgrid(1:size(ymesh2,2), 1:size(xmesh2,1), zmesh2);

    % gridded=10*log10(SI(xmesh3, ymesh3, zmesh3));

    figure("Units","inches","Position",[0,0,6,6])
    hold on



    mask = abs(elevations - 1.08)<0.1;
    [azmesh, rangemesh] = meshgrid(azimuths(mask), ranges);

    xmesh = rangemesh.*sind(azmesh);
    ymesh = rangemesh.*cosd(azmesh);
    % zmesh = rangemesh.*tand(elevations(mask));
    Zhhm = ZhhSur(:,mask);
    % Zhhm(mask2,:) = NaN;
    % Zhhm(Zhhm<0) = NaN;


    Zhhm(Zhhm<10) = NaN;
    figure("Units","inches","Position",[0,0,2.9,3])
    plotEarth(gca)
    hold on
    pcolor(xmesh/111000 + 120.434, ymesh/111000 + 22.527, Zhhm)
    shading flat
    xlim([min(xmesh/111000 + 120.434, [], 'all'), max(xmesh/111000 + 120.434, [], 'all')])
    ylim([min(ymesh/111000 +  22.527, [], 'all'), max(ymesh/111000 +  22.527, [], 'all')])
    yticks(21:0.5:24)
    xticks(119:0.5:122)
    yticklabels(string(yticks)+"^\circ N")
    xticklabels(string(xticks)+"^\circ E")
    x = cosd(0:360)
    y = sind(0:360)
    for rr = 50:50:150
        plot(x*rr/111 + 120.434, y*rr/111 + 22.527, 'k')
    end
    setMeteoColormap(gca, 'Zhh')
    daspect([1,1,1])
    print2(gcf, '/h/eol/nbarron/figures/publications/exmiras/TIMREX_Obs_SUR_Zhh.pdf')
    print2
    

    % pcolor(xmesh2(:,:,1), ymesh2(:,:,1), gridded(:,:,41))
    shading flat
    axis equal
    print2


    % % mask = abs(elevations - 1.08)<0.1 & abs(azimuths - 315)<45;
    
    

    % % nuimagesc(gca,ymesh, xmesh,  Zhhm)
    % setMeteoColormap(gca, 'Zhh')
    % ylim([-150000,150000])
    % xlim([-150000,150000])
    % shading flat
    % axis equal
    % print2
    
%}


%% interpolants/interpolate for the reflectivity and velocity fields
    Zhhsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), 10.^(double(Zhhm(:))./10), 'linear');
    Zdrsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(Zdrm(:)), 'linear');
    vrsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(VRm(:)), 'linear');
    Kdpsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(Kdpm(:)), 'linear');

    % interpolate to the model grid
    ZhhProfile = (10*log10(Zhhsi(rzmesh, zrmesh)));
    ZdrProfile = Zdrsi(rzmesh, zrmesh);
    vrProfile = vrsi(rzmesh, zrmesh);
    KdpProfile = Kdpsi(rzmesh, zrmesh);

    indMax = find(xGrid./meanVel > 4600+3600, 1,"first");

    ZhhProfile = ZhhProfile(1:indMax,:);
    ZdrProfile = ZdrProfile(1:indMax,:);
    vrProfile = vrProfile(1:indMax,:);
    KdpProfile = KdpProfile(1:indMax,:);

    mask =  ZdrProfile < -1 | ZhhProfile < 10;

    ZhhProfile(mask) = NaN;
    ZdrProfile(mask) = NaN;
    vrProfile(mask) = NaN;
    KdpProfile(mask) = NaN;

    ZhhProfile(:,zGrid<1000) = NaN;
    ZdrProfile(:,zGrid<1000) = NaN;
    vrProfile(:,zGrid<1000) = NaN;
%% perform the DSD assimilation
    mu = [];
    lambda = [];
    N0 = [];
    error = [];
    ZhhOpt = [];
    ZdrOpt = [];
    DmOpt = []
    inds = 2:size(ZhhProfile, 1)-1;
    ra = radar('S');
    figure
    hold on
    N = NaN(indMax-2, 41,250);
    parfor jj = 120:indMax-2
        % jj = 50
        jj
        indToTest = max(inds(jj)-1,1): min(inds(jj)+1, size(ZhhProfile,1)-2);
        
        [N0Temp, muTemp, lambdaTemp,fv,NTemp]=de.profileOptimizer(...
            ZhhProfile(indToTest,:)', ... 
            ZdrProfile(indToTest,:)', ...
            'KdpProfileObs',KdpProfile(indToTest,:)' ...
        );
        %% no Kdp method
        % [N0Temp]
        % [muTemp, lambdaTemp]
        % [N0Temp, muTemp, lambdaTemp,fv,NTemp]=de.profileOptimizer(...
        %     ZhhProfile(indToTest,:)', ... 
        %     ZdrProfile(indToTest,:)' ...
        %      ...
        % );
        % [N0Temp]
        % [muTemp, lambdaTemp]
        
    
        % subplot(1,3,1)
        % hold off
        % plot(ZhhProfile(indToTest,:)', zGrid)
        % hold on
        % [~, Zhhp, Zdrp] = de.estimateSingleRadarProfile(N0Temp, muTemp, lambdaTemp);
        % plot(Zhhp', zGrid)
        % subplot(1,3,2)
        % hold off
        % plot(ZdrProfile(indToTest,:)', zGrid)
        % hold on
        % plot(Zdrp', zGrid)
        % subplot(1,3,3)
        % hold off
        % plot(KdpProfile(indToTest,:)', zGrid)
        % hold on
        % [~,~, ~, Kdpp] = de.estimateSingleRadarProfile(N0Temp, muTemp, lambdaTemp);
        % plot(Kdpp', zGrid)
        % print2

        % figure
        % hold on
        % plot(de.D, squeeze(NTemp(end,:))')
        % print2

        % input('press enter to continue')
        if numel(NTemp) == 1
            NTemp = NaN(41,250);
        end
        N(jj,:,:) = NTemp;
        mu(jj) = muTemp;
        lambda(jj) = lambdaTemp;
        error(jj) = fv;
        N0(jj) = N0Temp;
        ZdrOpt(jj) = mean(ZdrProfile(indToTest,:), 'all', 'omitnan');
        ZhhOpt(jj) = mean(ZhhProfile(indToTest,:), 'all', 'omitnan');

        [ra.calcZhh(de.getNFromN0MuLambda(N0Temp, muTemp, lambdaTemp)), mean(ZhhProfile(indToTest,:), 'all', 'omitnan')]
        % fv
    end

    figure
    scatter(lambda(N0>1e5), mu(N0>1e5))
    hold on
    lambdai = linspace(0,20,100);
    mui = -0.016*lambdai.^2 + 1.213*lambdai - 1.957;% initial guess for mu
    plot(lambdai, mui, 'k--')
    mui = -0.0201*lambdai.^2 + 0.902*lambdai - 1.78 ;% initial guess for mu
    plot(lambdai, mui, 'r--')
    print2

%% Run EXMIRAS with the optimized DSD
    %% Set up the observation profiles
    obsRef = ZhhOpt;
    obsZdr = ZdrOpt;
    obsDm = DmOpt;
    obsSeconds = xGrid(1:indMax)./meanVel;

    % initialize exmiras object.
    ex = exmiras;

    % set up model domain
    ex.xgrid = [50];
    ex.ygrid = [50];
    ex.zgrid = 25:50:2025;

    % set starting state variables
    es = @(T) 2.53e9*exp(-5420./(T)); %hPa
    ex.T = TProfile;
    ex.p = presProfile;
    ex.pv = es(ex.T).*RHProfile;

    ex.u = zeros(size(ex.T));
    ex.v = zeros(size(ex.T));
    ex.w = zeros(size(ex.T));

    % reshape to match model grid
    ex.T = reshape(ex.T, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);
    ex.p = reshape(ex.p, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);
    ex.pv = reshape(ex.pv, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);
    ex.u = reshape(ex.u, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);
    ex.v = reshape(ex.v, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);
    ex.w = reshape(ex.w, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);

    % initialize the radar simulator --- this is required for both the radar and the dsdAssimilation class
    ex = ex.initializeRadarSimulator('S');

    % set the initial DSD based on the first observation
    iObs = 1;
    ex.N = zeros(numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid), numel(ex.D));
    ex.N(1,1,end,:) = ex.da.getNFromN0MuLambda(N0(iObs), mu(iObs), lambda(iObs));
    ex.NP = ex.N;

    %set time step/total integration time
    ex.dt = 0.5;
    ex.st = max(obsSeconds)./ex.dt;

    %% initialize the state/analysis variables
    sz = numel(ex.zgrid);
    numSteps = 6500*2;
    % numSteps = 6500*2;
    T = NaN(numSteps, sz);
    p = NaN(numSteps, sz);
    qv = NaN(numSteps, sz);
    pv = NaN(numSteps, sz);
    dBZhh = NaN(numSteps, sz);
    dBZvv = NaN(numSteps, sz);
    theta = NaN(numSteps, sz);
    RR = NaN(numSteps, sz);
    rhohv = NaN(numSteps, sz);
    kdp = NaN(numSteps, sz);
    % save initial temperature profile for later plotting
    T0 = ex.T;

    %% integrate forward in time, updating the model state with observations as they become available
    ex = ex.integrate;
    dts = ex.dt;
    N(isnan(N)) = 0;
    fig = figure('units','inches','position',[0,0,5.75,3]);
    hold on
    for i = 2:numSteps
        % check if we need to update the model with a new observation
        if (i-1)*ex.dt>=obsSeconds(iObs+1) & i*ex.dt<obsSeconds(iObs+2)
            % update the model with the next observation
            fprintf('updating model with observation at t = %4.0f seconds\n', obsSeconds(iObs+1))
            NTemp0 = ex.N;

            iObs = iObs + 1;
            % dbz = zeros(numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid));
            % dbz(end) = obsRef(iObs);
            % zdr = zeros(numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid));
            % zdr(end) = max(min(obsZdr(iObs),3),0.2);

            %! get the new DSD and add it to the top of the domain
            %! old method, may revert to this if problems arise
            %! Ntemp = ex.da.getNFromN0MuLambda(N0(iObs), mu(iObs), lambda(iObs));
            %! Ntemp(:,:,1:end-1,:) = 0;
            %! ex.NP = Ntemp;


            Ntemp = zeros(size(ex.N));
            for kk = 1:40
                % Ntemp(1,1,kk,:) = mean([reshape(ex.N(1,1,kk,:),[],1), reshape(N(iObs,kk,:),[],1)],2);
                Ntemp(1,1,kk,:) = N(iObs,kk,:);
            end
            Ntemp(1,1,end,:) = N(iObs,end,:);
            ex.N = Ntemp;
            Ntemp(:,:,1:end-1,:) = 0;
            ex.NP = Ntemp;
    
        end


        % integrate forward one time step
        ex = ex.integrate;

        % save the radar/state variables
        T(i,:) = squeeze(ex.T(1,1,:))-squeeze(T0(1,1,:));
        theta(i,:) = ex.theta(1,1,:);
        p(i,:) = ex.p(1,1,:);
        qv(i,:) = ex.qv(1,1,:);
        pv(i,:) = ex.pv(1,1,:);
        dBZhh(i,:) = ex.Zhh(1,1,:);
        dBZvv(i,:) = ex.Zvv(1,1,:);
        RR(i,:) = ex.RR(1,1,:);
        rhohv(i,:) = ex.rhohv(1,1,:);
        kdp(i,:) = ex.kdp(1,1,:);
        dts(i) = sum(dts) + ex.dt;
    end
    Zdr = dBZhh - dBZvv;

    save('./TIMREX-SPOL-EXMIRAS-withKdp.mat', 'T', 'p', 'qv', 'pv', 'dBZhh', 'dBZvv', 'Zdr', 'RR', 'rhohv', 'kdp', 'ex');

    load('./TIMREX-SPOL-EXMIRAS-withKdp.mat')
    


    %% Make plots of the results
    % load(sprintf('./SL-%d_%d.mat', levelOfAOSPRE,scanID), 'T', 'p', 'qv', 'pv', 'dBZhh', 'dBZvv', 'Zdr', 'RR', 'rhohv', 'kdp', 'ex');
    % [TT, ZT] = meshgrid((1:numSteps)*ex.dt, ex.zgrid);

    [~,dTdt]=gradient(T, ex.dt);
    dTdt = movmean(dTdt*3600,10*60,1);

    % TT = TT';
    % ZT = ZT';

    variablesToPlot = {dTdt, qv*100, dBZhh, Zdr, RR, rhohv, kdp};
    titles = {'dTdt', 'RH', 'Zhh', 'Zdr', 'RainRate', 'CC', 'Kdp'};
    
    numSteps = ex.st+60/ex.dt;
    [~,dT] = gradient(T, 0.5);
    dT = dT*3600; % convert to K/hr

    variablesToPlot = {dlog10(dT), 100*qv, dBZhh, Zdr, RR, rhohv, kdp};
    titles = {'dTdt', 'RH', 'Zhh', 'Zdr', 'RainRate', 'CC', 'KDP'};

    plotTimeProfile( ex.zgrid,(1:size(Zdr,1))*ex.dt,['/h/eol/nbarron/figures/publications/exmiras/TIMREX_test-'], variablesToPlot, titles)

    