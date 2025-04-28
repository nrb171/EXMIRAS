ev = evap;
ev.xgrid = [50];
ev.ygrid = [50];
ev.zgrid = 50:50:2050;
tab = readtable('/h/eol/nbarron/export/soundings/2015-06-14-DOC.txt');
tab.Properties.VariableNames = {'PRES', 'HGHT', 'TEMP', 'DWPT', 'RELH', 'MIXR', 'DRCT', 'SKNT', 'THTA', 'THTE', 'THTV'};

% temp = @(z) -9.8*z/1000 +  7.3402+273.15+tab.TEMP(1);
% pres = @(z) 1000*exp(-(z-51)/1000/8.4);
% es = @(T) 2.53e9*exp(-5420./(T)); %hPa

% dunion sounding
temp2 = [27.6, 26.7,21.6,17.5,9.4];
z2 = [0,138,823,1553,3190];
pres2 = [1016.5, 1000, 925, 850, 700];
RH = [77.2,78.9,81.7,69.5,33.5]
es2 = es(temp2+273.15).*RH/100;

z = 50:50:2050;
T = interp1(z2, temp2+273.15, ev.zgrid);
P = interp1(z2, pres2, ev.zgrid);
pv = interp1(z2, es2, ev.zgrid);


% ps = es(T);
% pv = [ps(z<tab.HGHT(1))*tab.RELH(1)/100, interp1(tab.HGHT, es(tab.TEMP+273.15), z(z>tab.HGHT(1)))/100];

% DWPT = tab.DWPT(1);
% DWPT(z>tab.HGHT(1)) = interp1(tab.HGHT, tab.DWPT, z(z>tab.HGHT(1)));


% skewt(tab.PRES,tab.TEMP, tab.DWPT);
% title('2015-06-14 00:00 UTC');
% print2(gcf, '~/figures/evap/sounding201506140000.png')



sx = numel(ev.xgrid);
sy = numel(ev.ygrid);
sz = numel(ev.zgrid);

ev.T = zeros(sx, sy, sz)
ev.T(1,1,:) = T;
Ti = T;
ev.p = zeros(sx, sy, sz);
ev.p(1,1,:) = P;
ev.pv = zeros(sx, sy, sz);
ev.pv(1,1,:) = pv;

ev.u = zeros(size(ev.T));
ev.v = zeros(size(ev.T));
ev.w = zeros(size(ev.T));


dBZi = 35 + rand(sx, sy, sz)*0;
Zdri = 0.8 + rand(sx, sy, sz)*0;
dBZi(:,:,1:size(dBZi,3)-1) = -inf;
Zdri(:,:,1:size(dBZi,3)-1) = 0;

ev = ev.initFromReflectivity(dBZi, Zdri);
ev.NP = ev.N;




ev.dt = 1;
ev.st = 600./ev.dt;
bandName = "S";





%% storm time?
numSteps = ev.st+600/ev.dt;
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
for i = 1:numSteps
    % tic
    ev = ev.integrate;
    T(i,:) = squeeze(ev.T(1,1,:))-squeeze(Ti)';
    theta(i,:) = ev.theta(1,1,:);
    p(i,:) = ev.p(1,1,:);
    qv(i,:) = ev.qv(1,1,:);
    pv(i,:) = ev.pv(1,1,:);
    dBZhh(i,:) = ev.Zhh(1,1,:);
    dBZvv(i,:) = ev.Zvv(1,1,:);
    RR(i,:) = ev.RR(1,1,:);
    rhohv(i,:) = ev.rhohv(1,1,:);
    kdp(i,:) = ev.kdp(1,1,:);
    % toc
    % r = 0.9951 + 0.02510*ev.D - 0.03644*ev.D.^2 + 0.005303*ev.D.^3 - 0.0002492*ev.D.^4;
end
Zdr = dBZhh - dBZvv;




[TT, ZT] = meshgrid((1:numSteps)*ev.dt, ev.zgrid);



variablesToPlot = {T, qv, dBZhh, Zdr, RR, rhohv, kdp};
titles = {'temperature perturbation', 'water vapor saturation', 'Z_{HH}', 'Z_{DR}', 'Rain Rate', 'Correlation Coefficient', 'Specific Differential Phase'};
cLabels = {'K', '', 'dBZ', 'dBZ', 'mm/hr', ' ', 'deg/km'};
clims = {[-8,8], [0.3,1], [-2*10/6,50], [-2,5], [0, 8], [0.99, 1], [0.5e-2, 1e-1]};
colormaps = {colortables('newblue'), cbrewer('seq','PuBu',64), colortables('RADAR32'), colortables('hyperspectral'), cbrewer('seq','RdPu',64), flipud(turbo(64)), (autumn(64))+[-0.5,0,0.5]};
for i = 1:numel(variablesToPlot)
    fig = figure("Units", "inches", "Position", [0,0,5.75,1.5]);

    pcolor(TT, ZT, real(variablesToPlot{i}'))
    ylabel('z [m]')
    % t = strsplit(titles{i}, '_');
    title(titles{i})
    cb = colorbar;
    cb.Label.String = cLabels{i};
    colormap(colormaps{i})
    clim(clims{i})
    shading flat
    titles{i} = strrep(strrep(titles{i}, "{", ''), "}", '');
    print2(fig, strrep(sprintf('~/figures/evap/%s_20150614_%s.png', bandName,titles{i}), ' ', '_'))
end

