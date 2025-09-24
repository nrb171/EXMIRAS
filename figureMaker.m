



%%% DSD estimator ideal%%%
    ZhhXlim = [38,42];
    ZdrXlim = [0.58, 0.68];

    %% Zhh profiles
    fig1 = figure('Units', 'inches', 'Position', [0,0,3,1.5]);
    ax1 = axes(fig1);
    grid(ax1, 'on')
    hold on

    %% Zdr profiles
    fig2 = figure('Units', 'inches', 'Position', [0,0,3,1.5]);
    ax2 = axes(fig2);
    grid(ax2, 'on')
    hold on

    %% n_{D_m\in[0.1,1.6]}(D)
    fig3 = figure('Units', 'inches', 'Position', [0,0,3,3]);
    ax3 = axes(fig3);
    grid(ax3, 'on')
    hold on

    fig4 = figure('Units', 'inches', 'Position', [0,0,3,1.5]);
    ax4 = axes(fig4);
    grid(ax4, 'on')
    hold(ax4, 'on')

    
    %% deta 


    zgrid = 25:50:2025;

    de = dsdEstimator;
    de.ZhhProfileObs = ones(1,41)*40;
    de.ZdrProfileObs = ones(1,41)*0.6;
    de.RHProfileObs = ones(1,41)*0.9;
    % de.RHProfileObs = RHProfile;
    de.bandName = 'C';
    % de = de.getNtop();
    
    % Zhh = [];
    % Zdr = [];
    % ZhhCalculated = [];
    colors = winter(numel(de.DmGrid));
    % colors2 = autumn(numel(de.DmGrid));

    D = logspace(log10(0.1), log10(8), 250);
    
    DmGrid = 0.1:0.3:1.6
    colors = winter(numel(DmGrid));
    meanD = [];
    maxD = 0.1:0.5:1.6;
    maxN = [];
    meanN = [];
    for jj = 1:numel(DmGrid)

        [dN, ZhhP, ZdrP, deta]=de.estimateSingleRadarProfile(de.ZhhProfileObs(end), de.ZdrProfileObs(end), DmGrid(jj));

        meanD(end+1) = trapz(D, dN(end,:).*D)./ trapz(D, dN(end,:));
        maxN(end+1) = max(dN(end,:));
        % maxN(end+1) = interp1(D, dN(end,:), maxD(end));
        meanN(end+1) = interp1(D, dN(end,:), meanD(end));

        deta = (dN(1,:)-dN(end,:))./max(dN(end,:));
        hold on
        

        % plot(ax1,ZhhP, zgrid, 'Color', colors(jj,:))
        % hold on
        % plot(ZhhCalculated(:,jj), zgrid, 'Color', colors2(jj,:))
        % subplot(1,2,2)
       
        plot(ax1,(ZhhP-(ZhhP(end)) + de.ZhhProfileObs(1)), zgrid, 'Color', colors(jj,:))
        % hold on
        % plot(ZhhCalculated(:,jj), zgrid, 'Color', colors2(jj,:))
        % subplot(1,2,2)
        plot(ax2,((ZdrP - (ZdrP(end))) + de.ZdrProfileObs(1)), zgrid, 'Color', colors(jj,:))

        plot(ax4, D, dlog10(squeeze(deta)*1000), 'Color', colors(jj,:))
    end
    grid(ax1, 'on')
    xlim(ax1, ZhhXlim)
    print2(fig1, '~/figures/publications/exmiras/dsdestimator-ideal-Zhh-vs-height-humid.pdf')

    grid(ax2, 'on')
    xlim(ax2, ZdrXlim)
    print2(fig2)
    print2(fig2, '~/figures/publications/exmiras/dsdestimator-ideal-Zdr-vs-height-humid.pdf')
    
    
    xlim(ax3, [0,4])
    % ylim(ax3, [1, 2e3])
    set(ax3, 'YScale', 'log')
    plot(ax3, maxD, maxN, 'k')
    % plot(ax3, meanD, meanN, 'r')
    print2(fig3, '~/figures/publications/exmiras/dsdestimator-ideal-N-vs-D.pdf')
    
    yticks(ax4, -3:1:3)
    ylim(ax4, [-3,3])
    yticklabels(ax4, ["-1", "-0.1", "-0.01", "0", "0.01", "0.1", "1"])
    xlim(ax4, [0,4])
    grid(ax4, 'on')
    % ylim(ax4, [-0.01, 0.01])
    print2(fig4, './.temp.png')
    % ylim(ax4, [-0.01, 0.01])
    print2(fig4, '~/figures/publications/exmiras/dsdestimator-ideal-deta-vs-D-humid.pdf')
    print2(fig2)

%% repeat, but for dry air

      %% Zhh profiles
    fig1 = figure('Units', 'inches', 'Position', [0,0,3,1.5]);
    ax1 = axes(fig1);
    grid(ax1, 'on')
    hold on

    %% Zdr profiles
    fig2 = figure('Units', 'inches', 'Position', [0,0,3,1.5]);
    ax2 = axes(fig2);
    grid(ax2, 'on')
    hold on

    %% n_{D_m\in[0.1,1.6]}(D)
    fig3 = figure('Units', 'inches', 'Position', [0,0,3,3]);
    ax3 = axes(fig3);
    grid(ax3, 'on')
    hold on

    fig4 = figure('Units', 'inches', 'Position', [0,0,3,1.5]);
    ax4 = axes(fig4);
    grid(ax4, 'on')
    hold(ax4, 'on')


    zgrid = 25:50:2025;

  
    de.RHProfileObs = ones(1,41)*0.3;
    % de.RHProfileObs = RHProfile;
    de.bandName = 'C';
    % de = de.getNtop();
    
    % Zhh = [];
    % Zdr = [];
    % ZhhCalculated = [];
    % colors = winter(numel(de.DmGrid));
    colors2 = autumn(numel(de.DmGrid));

    D = logspace(log10(0.1), log10(8), 250);
    
    DmGrid = 0.1:0.3:1.6
    colors = autumn(numel(DmGrid));
    meanD = [];
    maxD = 0.1:0.5:1.6;
    maxN = [];
    meanN = [];
    for jj = 1:numel(DmGrid)

        [dN, ZhhP, ZdrP, deta]=de.estimateSingleRadarProfile(de.ZhhProfileObs(end), de.ZdrProfileObs(end), DmGrid(jj));

        meanD(end+1) = trapz(D, dN(end,:).*D)./ trapz(D, dN(end,:));
        maxN(end+1) = max(dN(end,:));
        % maxN(end+1) = interp1(D, dN(end,:), maxD(end));
        meanN(end+1) = interp1(D, dN(end,:), meanD(end));

        deta = (dN(1,:)-dN(end,:))./max(dN(end,:));
        hold on
        

        plot(ax1,(ZhhP-(ZhhP(end)) + de.ZhhProfileObs(1)), zgrid, 'Color', colors(jj,:))
        % hold on
        % plot(ZhhCalculated(:,jj), zgrid, 'Color', colors2(jj,:))
        % subplot(1,2,2)
        plot(ax2,((ZdrP - (ZdrP(end))) + de.ZdrProfileObs(1)), zgrid, 'Color', colors(jj,:))

        plot(ax3, D, (squeeze(sum(dN(end,:),1))), 'Color', colors(jj,:))
        % plot(ax3, D, squeeze(sum(dN(1,:),1)), ':', 'Color', colors(jj,:))

        plot(ax4, D, dlog10(squeeze(deta)*1000), 'Color', colors(jj,:))
    end

    % xlim(ax2, [0.85,1])
    % xlim(ax1, [floor(min(ZhhP(:)))-0.5, ceil(max(ZhhP(:)))+0.5])
    xlim(ax1, ZhhXlim)
    grid(ax1, 'on')
    print2(fig1, '~/figures/publications/exmiras/dsdestimator-ideal-Zhh-vs-height-dry.pdf')
    grid(ax2, 'on')
    xlim(ax2, ZdrXlim)
    print2(fig2, '~/figures/publications/exmiras/dsdestimator-ideal-Zdr-vs-height-dry.pdf')
    
    
    % xlim(ax3, [0,4])
    % ylim(ax3, [1, 2e3])
    % set(ax3, 'YScale', 'log')
    % plot(ax3, maxD, maxN, 'k')
    % plot(ax3, meanD, meanN, 'r')
    % print2(fig3, '~/figures/publications/exmiras/dsdestimator-ideal-N-vs-D.pdf')
    
    yticks(ax4, -3:1:3)
    ylim(ax4, [-3,3])
    yticklabels(ax4, ["-1", "-0.1", "-0.01", "0", "0.01", "0.1", "1"])
    xlim(ax4, [0,4])
    grid(ax4, 'on')
    % ylim(ax4, [-0.01, 0.01])
    print2(fig4, './.temp.png')
    % ylim(ax4, [-0.01, 0.01])
    print2(fig4, '~/figures/publications/exmiras/dsdestimator-ideal-deta-vs-D-dry.pdf')
    print2(fig2)
%% DSD Estimator precip yonaguni 2022-06-16 12
    caseName = 'yonaguni-2022061612';
    ua=readmatrix("/h/eol/nbarron/workshop/EXMIRAS/test/data/yonaguni_2022061612.L2.csv");
    p = ua(:,3); % hPa
    T = ua(:,4); % C
    dwp = ua(:,12); % C
    RH = ua(:,5); % %
    Z = ua(:,10); % m


    zgrid = 25:50:2025;
    RHProfile = interp1(Z, RH, zgrid, "linear", "extrap")/100;

    de = dsdEstimator;
    de.ZhhProfileObs = ones(1,41)*50;
    de.ZdrProfileObs = ones(1,41)*1.5;
    % de.RHProfileObs = ones(1,41)*0.5;
    de.RHProfileObs = RHProfile;
    de.bandName = 'C';
    de = de.getNtop();
    de=de.estimateRadarProfiles();


    ex = exmiras;
    ex = ex.initFromLambdaName('C');
    ex.N = zeros(1,1,1,numel(ex.D));
    Zhh = [];
    Zdr = [];
    ZhhCalculated = [];
    colors = winter(numel(de.DmGrid));
    colors2 = autumn(numel(de.DmGrid));
    fig1 = figure('Units', 'inches', 'Position', [0,0,3,3]);
    ax1 = axes(fig1);
    hold on
    fig2 = figure('Units', 'inches', 'Position', [0,0,3,3]);
    ax2 = axes(fig2);
    hold on

    fig3 = figure('Units', 'inches', 'Position', [0,0,3,3]);
    ax3 = axes(fig3);
    hold on
    for jj = 1:numel(de.DmGrid)
        hold on
        

        plot(ax1,de.ZhhProfileSim(:,jj), zgrid, 'Color', colors(jj,:))
        % hold on
        % plot(ZhhCalculated(:,jj), zgrid, 'Color', colors2(jj,:))
        % subplot(1,2,2)
        plot(ax2,de.ZdrProfileSim(:,jj), zgrid, 'Color', colors(jj,:))

        plot(ax3, de.ex.D, squeeze(sum(de.dNProfile(:,jj,:),1)), 'Color', colors(jj,:))
    end

    xlim(ax2, [0,5])
    xlim(ax1, [40,55])
    print2(fig1, sprintf('~/figures/publications/exmiras/dsdestimator-%s-Zhh-vs-height.pdf', caseName))

    print2(fig2, sprintf('~/figures/publications/exmiras/dsdestimator-%s-Zdr-vs-height.pdf', caseName))


    xlim(ax3, [0,5])
    ylim(ax3, [-1e5, 1e5])
    print2(fig3, sprintf('~/figures/publications/exmiras/dsdestimator-%s-N-vs-D.pdf', caseName))


    figure('Units', 'inches', 'Position', [0,0,3,3])
    skewt(p, T, dwp)
    print2(sprintf('~/figures/publications/exmiras/dsdestimator-%s-skewt.pdf', caseName))



%% process/plot RHI for case
    d = "/h/eol/nbarron/work/precip/sea-pol/rhi/";
    caseName = 'yonaguni-2022061612';


    files = dir(d+"cfrad.20220616_154212.711_to_20220616_154417.869_SEAPOL_RHI.nc");


    % cft = CfRadialToolkit();
    % cft.plotCfRadial(files, ...
    %     'VariableName', ["DBZ"], ...
    %     'Units', 'inches', ...
    %     'Position', [0,0,3,1.5], ...
    %     'MeteoColormapName', "Zhh", ...
    %     'SaveDirectory', '~/figures/publications/exmiras/sea-pol-rhi/yonaguni-', ...
    %     'ylim', [0,12], ...
    %     'xlim', [0, 100], ...
    %     'Platform', 'SEA-POL', ...
    %     'Azimuth', [45:45:315]...
    % );
    % cft.plotCfRadial(files, ...
    %     'VariableName', ["ZDR"], ...
    %     'Units', 'inches', ...
    %     'Position', [0,0,3,1.5], ...
    %     'MeteoColormapName', "Zdr", ...
    %     'SaveDirectory', '~/figures/publications/exmiras/sea-pol-rhi/yonaguni-', ...
    %     'ylim', [0,12000], ...
    %     'xlim', [0, 100], ...
    %     'Platform', 'SEA-POL', ...
    %     'Azimuth', [45:45:315]...
    % );

    cft = CfRadialToolkit();
    range = ncread(fullfile(files(1).folder,files(1).name), 'range');
    azimuth = ncread(fullfile(files(1).folder,files(1).name), 'azimuth');
    elevation = ncread(fullfile(files(1).folder,files(1).name), 'elevation');
    Zhh = ncread(fullfile(files(1).folder,files(1).name), 'DBZ_TOT');
    % Zhh = ncread(fullfile(files(1).folder,files(1).name), 'DBZ_L2');
    Zdr = ncread(fullfile(files(1).folder,files(1).name), 'ZDR');
    vr = ncread(fullfile(files(1).folder,files(1).name), 'VEL');

    mask = abs(azimuth-135)<5; % only look at azimuths near 180 deg
    [elevmesh,rangemesh] = meshgrid(elevation(mask), range);
    xmesh = rangemesh.*cosd(elevmesh);
    ymesh = rangemesh.*sind(elevmesh);
    Zhhm = Zhh(:,mask);
    Zdrm = Zdr(:,mask);
    vrm = vr(:,mask);

    figure('Units', 'inches', 'Position', [0,0,3,1.5])
    pcolor(xmesh/1000, ymesh/1000, Zhhm)
    ylim([0,12])
    xlim([0,25])
    shading flat
    setMeteoColormap(gca, 'Zhh')
    print2(sprintf('~/figures/publications/exmiras/dsdestimator-%s-Zhh-rhi.pdf', caseName))
    % print2


    figure('Units', 'inches', 'Position', [0,0,3,1.5])
    pcolor(xmesh/1000, ymesh/1000, Zdrm)
    ylim([0,12])
    xlim([0,25])
    shading flat
    setMeteoColormap(gca, 'Zdr')
    print2(sprintf('~/figures/publications/exmiras/dsdestimator-%s-Zdr-rhi.pdf', caseName))



    Zhhsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(Zhhm(:)), 'linear');
    Zdrsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(Zdrm(:)), 'linear');
    vrsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(vrm(:)), 'linear');

    %% estimate tilt of precipitation, draw cuts along the tilt
    theta = atand(2.200); %

    startingPoints = 5000:50:8000;
    for ii = 1:numel(startingPoints)

        sp = startingPoints(ii);
        xr = sp + cosd(theta).*[0:50:4000/sind(theta)];
        yr = sind(theta).*[0:50:4000/sind(theta)];

        ZhhProfileObs(ii,:) = Zhhsi(xr, yr);
        ZdrProfileObs(ii,:) = Zdrsi(xr, yr);
        vrProfileObs(ii,:) = vrsi(xr, yr);
    end
  


    colors = summer(size(ZhhProfileObs,1));
    fig1 = figure('Units', 'inches', 'Position', [0,0,3,3]);
    ax1 = axes(fig1);
    hold on
    fig2 = figure('Units', 'inches', 'Position', [0,0,3,3]);
    ax2 = axes(fig2);
    hold on
    for jj = 1:size(ZhhProfileObs,1)
        plot(ax1,ZhhProfileObs(jj,:), yr, 'Color', colors(jj,:))
        plot(ax2,ZdrProfileObs(jj,:), yr, 'Color', colors(jj,:))
    end
    print2(fig1, sprintf('~/figures/publications/exmiras/dsdestimator-%s-Zhh-profile-vs-height.pdf', caseName))
    print2(fig2, sprintf('~/figures/publications/exmiras/dsdestimator-%s-Zdr-profile-vs-height.pdf', caseName))

    %% calculate CFAD of Zhh and Zdr
    zdrgridcfad = 0:0.5:5;
    zhhgridcfad = 20:2:56;

    zhhcfad = zeros(size(ZhhProfileObs,2), numel(zhhgridcfad)-1);
    zdrcfad = zeros(size(ZhhProfileObs,2), numel(zdrgridcfad)-1);
    for kk = 1:size(ZhhProfileObs,2)
        for ii = 1:numel(zhhgridcfad)-1
            zhhcfad(kk,ii) = sum(ZhhProfileObs(:,kk) >= zhhgridcfad(ii) & ZhhProfileObs(:,kk) < zhhgridcfad(ii+1));
        end
        for ii = 1:numel(zdrgridcfad)-1
            zdrcfad(kk,ii) = sum(ZdrProfileObs(:,kk) >= zdrgridcfad(ii) & ZdrProfileObs(:,kk) < zdrgridcfad(ii+1));
        end
    end

   

%% calculate difference from simulated to observed
    % load the observed profiles
    caseName = 'yonaguni-2022061612';
    ua=readmatrix("/h/eol/nbarron/workshop/EXMIRAS/test/data/yonaguni_2022061612.L2.csv");
    p = ua(:,3); % hPa
    T = ua(:,4); % C
    dwp = ua(:,12); % C
    RH = ua(:,5); % %
    Z = ua(:,10); % m


    zgrid = 25:50:2025;
    RHProfile = interp1(Z, RH, zgrid, "linear", "extrap")/100;

    de = dsdEstimator;
    de.bandName = 'C';
    de.RHProfileObs = RHProfile;
    % de.RHProfileObs = ones(size(RHProfile))*0.3;
    % ZhhProfileObs(:,yr >1500) = NaN;
    % ZdrProfileObs(:, yr >1500) = NaN;
    mask = 1:21;
    [ZhhOpt, ZdrOpt, DmOpt]=de.profileOptimizer(...
        interp1(yr, ZhhProfileObs(mask,:)', 25:50:2025), ... 
        interp1(yr, ZdrProfileObs(mask,:)', 25:50:2025) ...
    );


    mu = [];
    gamma = [];
    parfor ii = 1:size(ZhhProfileObs,1)-1

        mask = ii:ii+1;
        [ZhhOpt, ZdrOpt, DmOpt]=de.profileOptimizer(...
            interp1(yr, ZhhProfileObs(mask,:)', 25:50:2025), ... 
            interp1(yr, ZdrProfileObs(mask,:)', 25:50:2025) ...
        );
        ra = radar('C');
        [~, ~,muTemp, gammaTemp] = ra.initFromDm(ZhhOpt, ZdrOpt, DmOpt);
        mu(ii) = muTemp;
        gamma(ii) = gammaTemp;
    end

    figure('Units', 'inches', 'Position', [0,0,3,3])
    hold on
    scatter(gamma, mu)
    muFun = @(gamma) -0.016*gamma.^2 + 1.213*gamma - 1.957;
    gammaRange = linspace(0, 20, 100);
    plot(gammaRange, muFun(gammaRange), 'k--')

    muFun2 = @(gamma)-0.0201.*gamma.^2 + 0.902.*gamma -  1.718;
    plot(gammaRange, muFun2(gammaRange), 'r--')
    xlim([0,20])
    ylim([-2,14])
    print2
    
    

    fig3 = figure('Units', 'inches', 'Position', [0,0,3,3]);
    ax3 = axes(fig3);
    hold on
    colormap(jet(32))
    contourf((zhhgridcfad(1:end-1)+zhhgridcfad(2:end))/2, yr, zhhcfad./max(zhhcfad,[],2),32, 'LineStyle', 'none')
    ylim([0,2500])
    colors = winter(numel(0.3:0.3:1.6));
    print2(fig3, sprintf('~/figures/publications/exmiras/dsdestimator-%s-Zhh-cfad.pdf', caseName))
    for dd = 0.3:0.3:1.6
        [~,ZhhOptP,ZdrOptP] = de.estimateSingleRadarProfile(ZhhOpt, ZdrOpt, dd);
        plot(ax3,ZhhOptP, 25:50:2025, 'Color', colors(round((dd-0.1)/0.3)+1,:))
        % plot(ZhhOptP, 25:50:2025, 'Color', colors(round((dd-0.1)/0.3)+1,:))
    end
    for ii = mask
        plot(ax3, ZhhProfileObs(ii,:), yr)
    end
    xlim(ax3, [20, 56])
    [~,ZhhOptP,ZdrOptP] = de.estimateSingleRadarProfile(ZhhOpt, ZdrOpt, DmOpt);
    plot(ax3,ZhhOptP, 25:50:2025, 'k')
    print2(fig3, sprintf('~/figures/publications/exmiras/dsdestimator-%s-Zhh-cfad+optimizedFit.pdf', caseName))
    print2
    
    


    fig4 = figure('Units', 'inches', 'Position', [0,0,3,3]);
    ax4 = axes(fig4);
    hold on
    colormap(jet(32))
    contourf((zdrgridcfad(1:end-1)+zdrgridcfad(2:end))/2, yr, zdrcfad./max(zdrcfad,[],2),32, 'LineStyle', 'none')
    ylim([0,2500])
    print2(fig4, sprintf('~/figures/publications/exmiras/dsdestimator-%s-Zdr-cfad.pdf', caseName))

    colors = winter(numel(0.1:0.3:1.6));
    for dd = 0.1:0.3:1.6
        [~,ZhhOptP,ZdrOptP] = de.estimateSingleRadarProfile(ZhhOpt, ZdrOpt, dd);
        plot(ax4,ZdrOptP, 25:50:2025, 'Color', colors(round((dd-0.1)/0.3 )+1,:))
    end
    [N,ZhhOptP,ZdrOptP] = de.estimateSingleRadarProfile(ZhhOpt, ZdrOpt, DmOpt);
    plot(ax4,ZdrOptP, 25:50:2025, 'k')  
    print2
    print2(fig4, sprintf('~/figures/publications/exmiras/dsdestimator-%s-Zdr-cfad+optimizedFit.pdf', caseName))


   
    
%% plot CFADs for Zhh and Zdr


        
    d = "/h/eol/nbarron/work/precip/disdro/";

    files = dir(d+"*.nc");

    figure
    hold on
    for ii = 1:numel(files)
        f = files(ii).name;
        d = files(ii).folder;
        
        RR = ncread(fullfile(d,f), 'precip_rate');
        tt = ncread(fullfile(d,f), 'time'); 
        yyyy = str2double(f(14:17));
        mm = str2double(f(18:19));
        dd = str2double(f(20:21));
        dt = minutes(tt) + datetime(yyyy,mm,dd);
        % title(f)
        plot(dt, RR, '-k');
        % pause
    end

    xlim([datetime(2022, 6,16,12,0,0), datetime(2022, 6,16, 16,0,0)])
    print2

    d = "/h/eol/nbarron/work/precip/disdro/";
    f = "DSD_Yonaguni_20220616_qc.nc";

    RR = ncread(fullfile(d,f), 'precip_rate');
    [~,im]=max(RR);
    D = ncread(fullfile(d,f), 'particle_size');
    ndp = ncread(fullfile(d,f), 'qc_number_detected_particles');
    rs = ncread(fullfile(d,f), 'raw_spectrum');

    tt = ncread(fullfile(d,f), 'time');

    averageDropletSizeObs = trapz(D, sum(rs(:,:,im),2).*D)./trapz(D, sum(rs(:,:,im),2))
    % sum(rs(:,:,im),2))

    ndistro = zeros(size(D));
    for ii = -30:30
        ndistroTemp = sum(rs(:,:,im+ii),2);  


        ndistro = ndistro + ndistroTemp;
    end
    ndistronorm = ndistro./max(ndistro);

    figure('Units', 'inches', 'Position', [0,0,3,3])
    D2 = logspace(log10(0.1), log10(8), 250);
    hold on
    errors = [];
    colors = winter(numel(0.1:0.3:1.6));
    for dd = 0.4:0.3:1.3
        [N,ZhhOptP,ZdrOptP] = de.estimateSingleRadarProfile(ZhhOpt, ZdrOpt, dd);
        plot(D2, N(end,:)./max(N(end,:)), 'Color', colors(round((dd-0.1)/0.3 )+1,:))

        errors(end+1) = sum(abs(N(end,:)./max(N(end,:)) - interp1(D,ndistronorm,D2)));

        
    end

    [N,ZhhOptP,ZdrOptP] = de.estimateSingleRadarProfile(ZhhOpt, ZdrOpt, DmOpt);
    errors(end+1) = sum(abs(N(end,:)./max(N(end,:)) - interp1(D,ndistronorm,D2)));
    plot(D2, N(end,:)./max(N(end,:)),'k')
    % ylim([0, max(N(1,:))])
    % xlim([0 8])
    yyaxis right
    plot(D, ndistronorm)

    ylim([0, max(ndistronorm)])
    xlim([0.33 4])
    print2(sprintf('~/figures/publications/exmiras/dsdestimator-%s-N-vs-D+disdro.pdf', caseName))
    print2




%% plot time-height plots for TIMREX
    f = load('./TIMREX-SPOL-EXMIRAS.mat');
    ex = f.ex;
    T = f.T;
    p = f.p;
    qv = f.qv;
    dBZhh = f.dBZhh;
    dBZvv = f.dBZvv;
    Zdr = f.Zdr;
    RR = f.RR;
    rhohv = f.rhohv;
    kdp = f.kdp;

    % load('./KOKCcasestudy.mat');

    [TT, ZT] = meshgrid((1:numSteps)*ex.dt, ex.zgrid);

    [~,dTdt]=gradient(T, ex.dt);
    dTdt = movmean(dTdt*3600,10*60,1);
    TT = TT';
    ZT = ZT';
    variablesToPlot = {dTdt, qv*100, dBZhh, Zdr, RR, rhohv, kdp};
    titles = {'dTdt', 'RH', 'Zhh', 'Zdr', 'RainRate', 'CC', 'Kdp'};
    numSteps = ex.st+60/ex.dt;


    % scanID = 1;
    % load(sprintf('./SL-%d_%d.mat', levelOfAOSPRE,scanID), 'T', 'p', 'qv', 'pv', 'dBZhh', 'dBZvv', 'Zdr', 'RR', 'rhohv', 'kdp', 'ex');
    [~,dT] = gradient(T, 0.5);
    dT = dT*3600; % convert to K/hr
    variablesToPlot = {dlog10(dT), 100*qv, dBZhh, Zdr, RR, rhohv, kdp};
    titles = {'dTdt', 'RH', 'Zhh', 'Zdr', 'RainRate', 'CC', 'KDP'};
    plotTimeProfile( ex.zgrid,(1:numSteps)*ex.dt,['/h/eol/nbarron/figures/publications/exmiras/TIMREX_'], variablesToPlot, titles)

    averagedTdt = mean(dTdt(2000:10000,:), 'all');

%% plot time-height for ideal
    f = load('/h/eol/nbarron/work/exmirasDataS_ideal_300_0.30_40_0.50.mat');
    f = f.ExmirasRun;
    ex = f.ex;
    T = f.T;
    p = f.p;
    qv = f.qv;
    dBZhh = f.Zhh;
    dBZvv = f.Zvv;
    Zdr = f.Zdr;
    RR = f.RR;
    rhohv = f.rhohv;
    kdp = f.kdp;

    numsteps = size(T,1);

    [TT, ZT] = meshgrid(f.Grids.tgrid, f.Grids.zgrid);

    % [~,dTdt]=gradient(T, f.dt);

    % [~,dTdt]=gradient(T, ex.dt);
    % dTdt = movmean(dTdt*3600,10*60,1);
    % TT = TT';
    % ZT = ZT';
    % variablesToPlot = {dTdt, qv*100, dBZhh, Zdr, RR, rhohv, kdp};
    % titles = {'dTdt', 'RH', 'Zhh', 'Zdr', 'RainRate', 'CC', 'Kdp'};
    % numSteps = ex.st+60/ex.dt;


    % scanID = 1;
    % load(sprintf('./SL-%d_%d.mat', levelOfAOSPRE,scanID), 'T', 'p', 'qv', 'pv', 'dBZhh', 'dBZvv', 'Zdr', 'RR', 'rhohv', 'kdp', 'ex');
    [~,dT] = gradient(T, 0.5);
    dT = dT*3600; % convert to K/hr
    variablesToPlot = {dlog10(dT), 100*qv, dBZhh, Zdr, RR, rhohv, kdp};
    titles = {'dTdt', 'RH', 'Zhh', 'Zdr', 'RainRate', 'CC', 'KDP'};
    plotTimeProfile( f.Grids.zgrid,f.Grids.tgrid,['/h/eol/nbarron/figures/publications/exmiras/IDEAL_40_0.5'], variablesToPlot, titles)


%% plot tim-height for AOSPRE-SL simulation
    f = load('./AOSPRE-SL-EXMIRAS.mat');
    ex = f.ex;
    T = f.T;
    p = f.p;
    qv = f.qv;
    dBZhh = f.dBZhh;
    dBZvv = f.dBZvv;
    Zdr = f.Zdr;
    RR = f.RR;
    rhohv = f.rhohv;
    kdp = f.kdp;


    [TT, ZT] = meshgrid((1:numSteps)*ex.dt, ex.zgrid);

    [~,dTdt]=gradient(T, ex.dt);
    dTdt = movmean(dTdt*3600,10*60,1);
    TT = TT';
    ZT = ZT';
    variablesToPlot = {dTdt, qv*100, dBZhh, Zdr, RR, rhohv, kdp};
    titles = {'dTdt', 'RH', 'Zhh', 'Zdr', 'RainRate', 'CC', 'Kdp'};
    numSteps = ex.st+60/ex.dt;


    % scanID = 1;
    % load(sprintf('./SL-%d_%d.mat', levelOfAOSPRE,scanID), 'T', 'p', 'qv', 'pv', 'dBZhh', 'dBZvv', 'Zdr', 'RR', 'rhohv', 'kdp', 'ex');
    [~,dT] = gradient(T, 0.5);
    dT = dT*3600; % convert to K/hr
    variablesToPlot = {dlog10(dT), 100*qv, dBZhh, Zdr, RR, rhohv, kdp};
    titles = {'dTdt', 'RH', 'Zhh', 'Zdr', 'RainRate', 'CC', 'KDP'};
    plotTimeProfile( ex.zgrid,(1:numSteps)*ex.dt,['/h/eol/nbarron/figures/publications/exmiras/AOSPRE_SL'], variablesToPlot, titles)

    averagedTdt = mean(dTdt(1:7200,:), 'all', 'omitmissing');


%% plot lambda, mu vs. Dm contour plots
    da = dsdAssimilation('S');
    ra = radar('S');

    
    muRange = -2:0.1:15;
    lambdaRange = 0:0.1:20;
    Dm = NaN(numel(muRange), numel(lambdaRange));
    Dmean = NaN(numel(muRange), numel(lambdaRange));
    Zhh = NaN(numel(muRange), numel(lambdaRange));
    Zdr = NaN(numel(muRange), numel(lambdaRange));
    N02 = NaN(numel(muRange), numel(lambdaRange));
    for ii = 1:numel(muRange)
        for jj = 1:numel(lambdaRange)
            N0 = da.getN0FromZhhMuLambda(30, muRange(ii), lambdaRange(jj));
            N = da.getNFromN0MuLambda(N0, muRange(ii), lambdaRange(jj));
            [~, ind] = max(N);
            DmeanTemp = trapz(da.D, N.*da.D)./trapz(da.D, N);
            if DmeanTemp < 0.225 || DmeanTemp > 2
                continue
            end
            Dm(ii,jj) = da.D(ind);
            Dmean(ii,jj) = DmeanTemp;
            Zhh(ii,jj) = ra.calcZhh(N);
            Zdr(ii,jj) = ra.calcZdr(N);
            N02(ii,jj) = log10(N0);
        end
    end

    figure('Units', 'inches', 'Position', [0,0,1.9,2])
    ax = axes;
    hold on
    colormap(turbo)
    contourf(lambdaRange, muRange, Dmean, 32, 'LineStyle', 'none')
    muFun = @(gamma) -0.016*gamma.^2 + 1.213*gamma - 1.957;
    gammaRange = linspace(0, 20, 100);
    plot(gammaRange, muFun(gammaRange), 'k--')

    muFun2 = @(gamma)-0.0201.*gamma.^2 + 0.902.*gamma -  1.718;
    plot(gammaRange, muFun2(gammaRange), 'r--')
    % colorbar
    % xlabel('\Lambda (mm^{-1})')
    % ylabel('\mu')
    % title('\bar{D} (mm)')
    xlim([0,20])
    ylim([-2,15])
    
    [lambdaMesh, muMesh] = meshgrid(lambdaRange, muRange);
    pf=fit([muMesh(:), lambdaMesh(:)], Dmean(:), 'linearinterp')

    inds1 = find(abs(Dmean - 2) < 0.01)
    inds2 = find(abs(Dmean - 0.225) < 0.01)

    polytop = [2.3277, -0.9189];
    polybot = [0.2796, -2.0];

    print2('~/figures/publications/exmiras/lambda-mu-Dmean.pdf')

    figure('Units', 'inches', 'Position', [0,0,1,2])
    ax = axes;
    colormap(turbo)
    clim([0.225,2])
    colorbar
    ax.Visible = 'off';
    print2('~/figures/publications/exmiras/lambda-mu-Dmean-colorbar.pdf')


    figure('Units', 'inches', 'Position', [0,0,1.9,2])
    ax = axes;
    contourf(lambdaRange, muRange, N02, 32, 'LineStyle', 'none')
    colormap(winter)
    clim([-2,12])
    % colorbar
    % xlabel('\Lambda (mm^{-1})')
    % ylabel('\mu')
    % title('N_0 (log_{10}(m^{-3} mm^{-1}))')
    print2('~/figures/publications/exmiras/lambda-mu-N0.pdf')
    
    figure('Units', 'inches', 'Position', [0,0,1,2])
    ax = axes;
    colormap(winter)
    clim([-2,12])
    ax.Visible = 'off';
    cb=colorbar;
    cb.TickLabels = {'10^{-2}', '10^{0}', '10^{2}', '10^{4}', '10^{6}', '10^{8}', '10^{10}', '10^{12}'};
    print2('~/figures/publications/exmiras/lambda-mu-N0-colorbar.pdf')

    figure('Units', 'inches', 'Position', [0,0,1.9,2])
    ax = axes;
    contourf(lambdaRange, muRange, Zdr, 32, 'LineStyle', 'none')
    setMeteoColormap(gca, 'Zdr')
    % colorbar
    % xlabel('\Lambda (mm^{-1})')
    % ylabel('\mu')
    % title('Z_{dr} (dB)')
    print2('~/figures/publications/exmiras/lambda-mu-Zdr.pdf')


    
%% plot idealized simulation for the DSD estimation
    % load("~/work/exmirasData/LUTs-N0MuLambda-0.10_-1.76_0.00_1.00.mat")
    % da = dsdAssimilation('S');
    % ra = radar('S');
    % Zhh = NaN(size(ExmirasRun.N,1), size(ExmirasRun.N,2));
    % Zdr = NaN(size(ExmirasRun.N,1), size(ExmirasRun.N,2));
    % for ii = 1:size(ExmirasRun.N,1)
    %     for jj = 1:size(ExmirasRun.N,2)
            
    %         Zhh(ii,jj) = ra.calcZhh(ExmirasRun.N(ii,jj,:));
    %         Zdr(ii,jj) = ra.calcZdr(ExmirasRun.N(ii,jj,:));
    %     end
    % end

    % figure
    % pcolor(1:size(ExmirasRun.N,1), ExmirasRun.Grids.zgrid, Zhh')
    % colorbar
    % shading flat
    % setMeteoColormap(gca, 'Zhh')
    % print2

%% plot PPI/SUR scan for TiMREX case
    folder = "/h/eol/nbarron/work/timrex/spol/20080614/";
    filename = "cfrad.20080614_083002.000_to_20080614_083138.000_SPOLRVP8_v25_SUR.nc";
    % filename = "cfrad.20080614_091644.000_to_20080614_091820.000_SPOLRVP8_v45_SUR.nc";
    % filename = "cfrad.20080614_110733.000_to_20080614_111447.000_SPOLRVP8_v62_SUR.nc";

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

%% plot RHI for TiMREX case
        d = "/h/eol/nbarron/work/timrex/spol/20080614/";
        files = dir(d+"*RHI*");

        zGrid = 25:50:2025;
        xGrid = 0:150:150000;

        dz=100;
        [zrmesh, rzmesh] = meshgrid(zGrid, xGrid);

        Dmean = 1.1;
        vmean = -0.1021 + 4.932*Dmean - 0.955*Dmean.^2 + 0.07934*Dmean.^3 - 0.002362*Dmean.^4;

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


%% plot sounding for TIMREX
    zGrid = 25:50:2025;
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
    