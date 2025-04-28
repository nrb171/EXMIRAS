% KOKC case study
%% load radar temporal profiles
ref = ncread('./test/KTLXProfiles.nc', 'refProfile');
zdr = ncread('./test/KTLXProfiles.nc', 'zdrProfile');
% zdr(isnan(zdr)) = 0;
TTobs = ncread('./test/KTLXProfiles.nc', 'time');
TTobs = datetime(1970,1,1) + seconds(TTobs);
ZTobs = ncread('./test/KTLXProfiles.nc', 'height');
% plot profiles
variablesToPlot = {ref, zdr};
titles = {'Z_{HH}', 'Z_{DR}'};
cLabels = {'dBZ', 'dBZ'};
clims = {[-2*10/6,50], [-2,5]};
colormaps = {colortables('RADAR32'), colortables('hyperspectral')};
for i = 1:numel(variablesToPlot)
    fig = figure("Units", "inches", "Position", [0,0,5.75,1.5]);

    pcolor(TTobs, ZTobs, real(variablesToPlot{i}))
    ylabel('z [m]')
    % t = strsplit(titles{i}, '_');
    title(titles{i})
    cb = colorbar;
    cb.Label.String = cLabels{i};
    colormap(colormaps{i})
    clim(clims{i})
    shading flat
    xlim([datetime(2015,6,13,2,0,0), datetime(2015,6,13,6,0,0)])
    ylim([0,2500])
    titles{i} = strrep(strrep(titles{i}, "{", ''), "}", '');
    print2(fig, strrep(sprintf('./KOKC_%s.png', titles{i}), ' ', '_'))
end

% somewhat difficult to read, might be better to use a CFAD.
figure('units','inches','position',[0,0,5.75,1.5])
plot(TTobs, movmean(ref(8,:),3))
ylabel('reflectivity [dBZ]')
grid on
yyaxis right
plot(TTobs, movmean(zdr(8,:),3))
ylabel('differential reflectivity [dBZ]')
grid on
% title(sprintf('KOKC reflectivity and differential reflectivity at %s km', num2str(ZT(8)/1000)))
xlim([datetime(2015,6,13,2,0,0), datetime(2015,6,13,6,0,0)])
print2(gcf,'./KOKCprofiles.png', 'quality', '-r300')

% set the time period for the model run (in indices)
idxStart = 35;
idxEnd = 35+3.5*12+1;

obsSeconds = seconds(TTobs(idxStart:idxEnd) - TTobs(idxStart));
obsRef = mean(ref(7:8,idxStart:idxEnd),1);
obsZdr = zdr(8,idxStart:idxEnd);
obsZdr(obsZdr<0) = 0.1;
%subsample the data to smooth
obsRef = real(10*log10(interp1(obsSeconds, 10.^(obsRef/10), 0:60:seconds(TTobs(idxEnd)-TTobs(idxStart)), 'pchip', 'extrap')));
obsZdr = real(10*log10(interp1(obsSeconds, 10.^(obsZdr/10), 0:60:seconds(TTobs(idxEnd)-TTobs(idxStart)), 'pchip', 'extrap')));
obsSeconds = 0:60:seconds(TTobs(idxEnd)-TTobs(idxStart));



%% load sounding data from KOUN
tab = readmatrix('./test/data/KOUN-2015-06-13_00:00.txt');

height = tab(:,15);
agl = height - height(1);

temp = tab(:,3);
pres = tab(:,2);
dewpt = tab(:,4);
relh = tab(:,5);


fig = figure('units','inches','position',[0,0,3.25,3]);
skewt(pres, temp, dewpt);
print2('./KOKCsounding.png')

%% set up model run
% initialize exmiras object.
ex = exmiras;

% set up model domain
ex.xgrid = [50];
ex.ygrid = [50];
ex.zgrid = 50:50:2500;

% set starting state variables
es = @(T) 2.53e9*exp(-5420./(T)); %hPa
ex.T = interp1(agl, temp+273.15, ex.zgrid);
ex.p = interp1(agl, pres, ex.zgrid);
ex.pv = interp1(agl, es(temp+273.15).*relh/100, ex.zgrid);

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


dBZi = obsRef(1) + rand(numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid))*0;
Zdri = obsZdr(1) + rand(numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid))*0;
dBZi(:,:,1:size(dBZi,3)-1) = -inf;
Zdri(:,:,1:size(dBZi,3)-1) = 0;

% initialize dsd from reflectivity/zdr
ex = ex.initFromReflectivity(dBZi, Zdri);
ex.NP = ex.N;

%set time step/total integration time
ex.dt = 1;
ex.st = 3.5*3600./ex.dt;

% set radar band
bandName = "S";


sz = numel(ex.zgrid);
%% storm time?
numSteps = ex.st+60/ex.dt;
T = zeros(numSteps, sz);
p = zeros(numSteps, sz);
qv = zeros(numSteps, sz);
pv = zeros(numSteps, sz);
dBZhh = zeros(numSteps, sz);
dBZvv = zeros(numSteps, sz);
theta = zeros(numSteps, sz);
RR = zeros(numSteps, sz);
rhohv = zeros(numSteps, sz);
kdp = zeros(numSteps, sz);

% save initial temperature profile for later plotting
T0 = ex.T;

ex = ex.integrate;
dts = [];
dts = TTobs(idxStart);
iObs = 1;
for i = 2:numSteps
    if (i-1)*ex.dt>=obsSeconds(iObs+1) & i*ex.dt<obsSeconds(iObs+2)
        % update the model with the next observation
        fprintf('updating model with observation')
        NTemp0 = ex.N;

        iObs = iObs + 1;
        dbz = zeros(numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid));
        dbz(end) = obsRef(iObs);
        zdr = zeros(numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid));
        zdr(end) = obsZdr(iObs);
        ex = ex.initFromReflectivity(dbz, zdr);
        ex.N(1,1,1:end-1,:) = NTemp0(1,1,1:end-1,:);
        Ntemp = ex.N;
        Ntemp(:,:,1:end-1,:) = 0;
        ex.NP = Ntemp;


        %% plot the model output
        

        % r = 0.9951 + 0.02510*ex.D - 0.03644*ex.D.^2 + 0.005303*ex.D.^3 - 0.0002492*ex.D.^4;
        fprintf('updated model with observation')
        iObs = iObs + 1;
   
    end


    ex = ex.integrate;
    
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
    dts(i) = dts(i-1) + seconds(ex.dt);
end
Zdr = dBZhh - dBZvv;



[TT, ZT] = meshgrid((1:numSteps)*ex.dt, ex.zgrid);
variablesToPlot = {T, qv, dBZhh, Zdr, RR, rhohv, kdp};
titles = {'temperature perturbation', 'water vapor saturation', 'Z_{HH}', 'Z_{DR}', 'Rain Rate', 'Correlation Coefficient', 'Specific Differential Phase'};
cLabels = {'K', '', 'dBZ', 'dBZ', 'mm/hr', ' ', 'deg/km'};
clims = {[-12,12], [0.3,1], [-2*10/6,50], [-2,5], [0, 12], [0.99, 1], [0.5e-2, 0.5]};
colormaps = {colortables('newblue'), cbrewer('seq','PuBu',64), colortables('RADAR32'), colortables('hyperspectral'), cbrewer('seq','RdPu',64), flipud(turbo(64)), flipud(autumn(64))+[0,0,0.8]};
plotTimeProfile(dts, ex.zgrid,'./KOKC', variablesToPlot, titles, cLabels, clims, colormaps)