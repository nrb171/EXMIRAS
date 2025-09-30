d = "/h/eol/nbarron/work/precip/sea-pol/rhi/";
files = dir(d+"*cfrad*");

zGrid = 25:50:2025;
xGrid = 3000:50:20000;

dz=100;
[zrmesh, rzmesh] = meshgrid(zGrid, xGrid);

Dmean = 1.1;
vmean = -0.1021 + 4.932*Dmean - 0.955*Dmean.^2 + 0.07934*Dmean.^3 - 0.002362*Dmean.^4;


%% build map of sounding times
Soundings = struct();
Soundings.files = dir("/h/eol/nbarron/work/precip/upper-air/yonaguni*.csv");
Soundings.timesInString = string({Soundings.files.name});
for ii = 1:numel(Soundings.timesInString)
    Soundings.timesInString(ii) = extractBetween(Soundings.timesInString(ii), "yonaguni_", ".L2.csv");
    Soundings.timesInDatetime(ii) = datetime(Soundings.timesInString(ii), 'InputFormat', 'yyyyMMddHH');
    ua=readmatrix(fullfile(Soundings.files(ii).folder, Soundings.files(ii).name), 'NumHeaderLines', 46);
    p = ua(:,3); % hPa
    T = ua(:,4); % C
    dwp = ua(:,12); % C
    RH = ua(:,5); % %
    Z = ua(:,10); % m

    % zgrid = 25:50:2025;
    [~,uidx]= unique(Z);
    RHInterpolant = griddedInterpolant(Z(uidx), RH(uidx)/100, 'linear', 'linear');
    Soundings.RHProfile(ii,:) = RHInterpolant(zGrid);

    Soundings.LCL(ii) = 125*(T(2)-dwp(2)); % approximate LCL in m
end

Soundings.interpolantDatenum = griddedInterpolant(datenum(Soundings.timesInDatetime), Soundings.RHProfile, 'linear', 'linear');
Soundings.interpolantDatetime = @(dt) Soundings.interpolantDatenum(datenum(dt));



figure('Units', 'inches', 'Position', [0,0,3,3])
gammaEdges = 0:0.25:20;
        muEdges = -2:0.25:15;
        N = zeros(numel(muEdges)-1, numel(gammaEdges)-1);
for ii = 1:numel(files)

    %% load necessary rhi data
    range = ncread(fullfile(files(ii).folder,files(ii).name), 'range');
    azimuth = ncread(fullfile(files(ii).folder,files(ii).name), 'azimuth');
    elevation = ncread(fullfile(files(ii).folder,files(ii).name), 'elevation');
    Zhh = ncread(fullfile(files(ii).folder,files(ii).name), 'DBZ_TOT');

    if max(Zhh(:)) < 30
        continue
    end

    Zdr = ncread(fullfile(files(ii).folder,files(ii).name), 'ZDR');
    vr = ncread(fullfile(files(ii).folder,files(ii).name), 'VEL');

    dt = datetime(files(ii).name(7:21), 'InputFormat', 'yyyyMMdd_HHmmss');

    RHProfile = Soundings.interpolantDatetime(dt);

    de = dsdEstimator;
    de.bandName = 'C';
    de.RHProfileObs = RHProfile;
    

    for az = unique(round(azimuth/5)*5)'
        mask = abs(azimuth-az)<5; % only look at azimuths near 180 deg
        [elevmesh,rangemesh] = meshgrid(elevation(mask), range);
        xmesh = rangemesh.*cosd(elevmesh);
        ymesh = rangemesh.*sind(elevmesh);
        Zhhm = Zhh(:,mask);
        Zdrm = Zdr(:,mask);
        vrm = vr(:,mask);

        Zhhsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(Zhhm(:)), 'linear');
        Zdrsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(Zdrm(:)), 'linear');
        vrsi = scatteredInterpolant(double(xmesh(:)), double(ymesh(:)), double(vrm(:)), 'linear');


        % calculate the displacement of an average drop based on Doppler velocity
        uz = mean(vrsi(rzmesh, zrmesh),1, 'omitnan');
        dt = dz/vmean;
        dxProfile = cumsum(uz.*dt, "omitmissing")/10; 

        %% reconstruct the profiles to match shear
        for rr = 1:size(rzmesh,1)

            xGridToSample = rzmesh(rr,1) - dxProfile;
            zGridToSample = zrmesh(rr,:);

            ZhhProfile(rr,:) = Zhhsi(xGridToSample, zGridToSample);
            ZdrProfile(rr,:) = Zdrsi(xGridToSample, zGridToSample);
        end
        
        
        %% perform the DSD estimation only when the Zhh is above 30 dBZ
        inds = find(ZhhProfile(:,end)>30);

        mu = [];
        gamma = [];
        parfor jj = 1:numel(inds)
            indToTest = max(inds(jj)-1,1): min(inds(jj)+1, size(ZhhProfile,1));
            
            [ZhhOpt, ZdrOpt, DmOpt]=de.profileOptimizer(...
                ZhhProfile(indToTest,:)', ... 
                ZdrProfile(indToTest,:)' ...
            );
            ra = radar('C');
            ra.rngToggle = true; % toggle for random number generator
            [~, ~,muTemp, gammaTemp] = ra.initFromDm(ZhhOpt, ZdrOpt, DmOpt);
            mu(jj) = muTemp;
            gamma(jj) = gammaTemp;
        end

        %% create mu/gamma database
        save(sprintf('/h/eol/nbarron/work/dsd-estimator/%s-%1.0f.mat', files(ii).name(7:21), az), 'mu', 'gamma', 'dt')

        
        N = N+histcounts2(gamma,mu,gammaEdges, muEdges)';
        
    end
end

colormap(cbrewer('seq', 'YlGnBu', 256))



files = dir('/h/eol/nbarron/work/dsd-estimator/*.mat');
mu = [];
gamma = [];
for ii = 1:numel(files)
    f = load(fullfile(files(ii).folder, files(ii).name));
    mu = [mu, f.mu];
    gamma = [gamma, f.gamma];
end

figure('Units', 'inches', 'Position', [0,0,3,3])
N = histcounts2(gamma, mu, gammaEdges, muEdges)';
[gammaMesh, muMesh] = meshgrid((gammaEdges(1:end-1) + gammaEdges(2:end))/2, (muEdges(1:end-1) + muEdges(2:end))/2);
gammaRescaled = [];
muRescaled = [];
for jj = 1:numel(gammaMesh)
    N2 = round(log10(1+N(jj)));
    if N2 == 0 
        continue
    end
    gammaRescaled(end+1:end+N2) = gammaMesh(jj);
    muRescaled(end+1:end+N2) = muMesh(jj);
end


figure('Units', 'inches', 'Position', [0,0,3,3])
hold off
pcolor((gammaEdges(1:end-1) + gammaEdges(2:end))/2, (muEdges(1:end-1) + muEdges(2:end))/2, log10(N), 'EdgeColor', 'none')
hold on
plot(gammaRange, muFun(gammaRange), 'k--')
muFun = @(gamma) -0.016*gamma.^2 + 1.213*gamma - 1.957;
gammaRange = linspace(0, 20, 100);

muFun2 = @(gamma)-0.0201.*gamma.^2 + 0.902.*gamma -  1.718;
plot(gammaRange, muFun2(gammaRange), 'r--')

gammaRescaled(muRescaled < 0) = [];
muRescaled(muRescaled < 0) = [];
muRescaled(gammaRescaled < 1) = [];
gammaRescaled(gammaRescaled < 1) = [];

pf = polyfit(gammaRescaled, muRescaled, 2)
errorfun = @(x) abs((mean(nonzeros(((polyval(x, gammaMesh) - muMesh).^2).*log10(1+N)))).^(1/2));
error = rms((polyval(pf, gammaMesh) - muMesh).*(N), "all")
pf=fminsearchbnd(errorfun, [0,0,0], [-0.25, -3, -3], [0.25, 3, 3])
% pf=fminsearchbnd(errorfun, [0,0,0])
plot(gammaRange, polyval(pf, gammaRange), 'Color',[94,79,162]/255, 'LineWidth', 1.5)
xlim([0,20])
ylim([-2,14])   
colormap(cbrewer('seq', 'YlGnBu', 256))

otherFit = fit(gammaRescaled', muRescaled', 'power2');
error = rms((otherFit(gammaMesh) - muMesh(:)).*(N(:)), "all")
plot(gammaRange, otherFit(gammaRange), 'Color',[0.9290, 0.6940, 0.1250], 'LineWidth', 1.5)
print2('~/figures/publications/exmiras/dsdEstimatorMuvGamma.pdf')
xlabel('$\Gamma$')
ylabel('$\mu$')
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
