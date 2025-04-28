a = readtable('KOKC201506.csv');
dts = datetime(string(table2array(a(:,'date_time')))) + hours(6);
temps = table2array(a(:,'air_temperature'));
dewpts = table2array(a(:,'dew_point_temperature'));
pressures = table2array(a(:,'air_pressure_at_sea_level'));

fig = figure("Units", "inches", "Position", [0,0,5.75,1.5]);
plot(dts, temps+273.15, 'k')
hold on 
plot(dts, dewpts+273.15, 'b')
ylabel('temperature [K]')

mask = contains(string(table2array(a(:,'current_wx1'))), ["RA", "TS"]);
bar(dts(mask), 1.1*max(temps(mask)+273.15)*ones([1,nnz(mask)]), 1, 'FaceColor', [0.5,0.5,0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5)


grid on
ax = gca;
ax.XMinorTick = 'on';
ax.XAxis.MinorTickValues = datetime(2015,6,1)+hours(0:1:23*30);
ylim([min(temps)+273.15 - 0.5, max(temps+273.15)+0.5])
yyaxis right
plot(dts, pressures)
ylabel('pressure [hPa]')
legend('temperature', 'dewpoint', 'rain', 'pressure')

xlim([datetime(2015,6,1), datetime(2015,6,31)])


print2(fig, '~/figures/evap/KOKC201506.png', 'quality','-r300')


%remove diurnal cycle
Fs = 1/(5*60)
T = 1/Fs;             % Sampling period       
L = length(temps)             % Length of signal
t = (0:L-1)*T;        % Time vector
temps = fillmissing(temps, 'spline');

% fList = Fs/L*(0:L-1)*24*3600
fListShifted = Fs/L*(-(L/2):(L/2-1))*24*3600;

Y = fft(temps-movmean(temps,60*12/5));
% Y2 = fftshift(Y);

figure
% plot(fList,abs(Y2))
plot(fListShifted,fftshift(real(Y)))
% yscale('log')
print2()



%reconstruct diurnal cycle
Y3 = zeros(size(Y));
mask = abs(fListShifted) == min(abs(fListShifted)) | abs(abs(fListShifted)-1) <= min(abs(abs(fListShifted)-1));
mask = ifftshift(mask);
[i,mask] = maxk(abs(Y),2);
Y3(mask) = Y(mask);

X = real(ifft(Y3));

figure
plot(dts, X)
hold on
plot(dts, movmean(temps-movmean(temps,60*12/5),3))
xlim([datetime(2015,6,11), datetime(2015,6,15)])
print2()


%%%%
temps2 = interp1(dts, temps, dts(1):hours(1):dts(end),'cubic')
dts2 = dts(1):hours(1):dts(end);

figure
plot(dts2, temps2)
print2()

Y = fft(temps2 - movmean(temps2, 48));
% Y2 = fftshift(Y);
Fs = 1;
L = length(temps2);
fList = Fs/L*(-L/2:L/2-1)*24;

figure
plot(fList,fftshift(abs(Y)))
print2()

%remove diurnal cycle
Y3 = zeros(size(Y));
[~,mask] = maxk(abs(Y),1);
mask = abs(ifftshift(fList))==1;
Y3(mask) = Y(mask);

X = real(ifft(Y3));

figure

mask = contains(string(table2array(a(:,'current_wx1'))), ["RA", "TS"]);
bar(dts(mask), 1.1*max(temps(mask)+273.15)*ones([1,nnz(mask)]), 1, 'FaceColor', [0.5,0.5,0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5)
hold on
plot(dts2, X+movmean(temps2, 48))
hold on
ylim([20,35])
xlim([datetime(2015,6,11), datetime(2015,6,15)])
plot(dts, temps)
plot(dts2, movmean(temps2, 48))




grid on
% ax = gca;
% ax.XMinorTick = 'on';
% ax.XAxis.MinorTickValues = datetime(2015,6,1)+hours(0:1:23*30);
ylim([min(temps) - 0.5, max(temps)+0.5])
yyaxis right
% plot(dts2, temps2-(X+movmean(temps2, 48)))
print2(gcf,'~/figures/evap/decomp.png')




figure('units','inches','position',[0,0,5.75,1.5])

mask = contains(string(table2array(a(:,'current_wx1'))), ["RA", "TS"]);
bar(dts(mask), 1.1*max(temps(mask)+273.15)*ones([1,nnz(mask)]), 1, 'FaceColor', [0.5,0.5,0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5)
hold on
plot(dts2, movmean(temps2, 48)+X, "LineWidth", 1.2,"Color", "k")
colors = gradient(-X-movmean(temps2, 48))+gradient(movmean(temps2,2));
colorsInterp = interp1(dts2, colors, dts);
cmap = flipud(cbrewer('div', 'Spectral', 64));
cmapInterpolant = griddedInterpolant(linspace(-1.5,1.5,64), cmap,"nearest");
colormap(cmap)
cb=colorbar
cb.Label.String = 'temperature rate difference [C hr^{-1}]';
% set(cb, 'Label', 'temperature rate difference [C hr^{-1}]')
clim([-1.5,1.5])
xlim([datetime(2015,6,13), datetime(2015,6,14)])
% 
temps3 = movmean(temps, 12);
% plot(dts, temps3, 'k', 'LineWidth', 2.5)

% p = plot(dts, temps3, 'LineWidth', 1);
% set(p.Edge, 'ColorBinding', 'interpolated', 'ColorData', uint8([cmapInterpolant(colorsInterp)*255, ones([numel(colorsInterp),1])*255]'))
X2 = movmean(temps2, 48);
for i = 1:numel(dts)-1
    % plot(dts2(i:i+1), [X(i)+X2(i), X(i+1)+X2(i+1)], 'Color', cmapInterpolant(colors(i)), 'LineWidth', 2)
    xlims = get(gca, 'XLim');
    if dts(i) < xlims(2) && dts(i+1) > xlims(1)
        plot(dts(i:i+1)+minutes([0,5]), temps3(i:i+1), 'Color', cmapInterpolant(mean(colorsInterp(i:i+1))), 'LineWidth', 1.5)
    end
    % hold on
end

plot(dts, movmean(dewpts,12), 'b')
ylim([18,35])
ylabel('temperature [C]')
yyaxis right
plot(dts, movmean(pressures,12))
ylabel('pressure [hPa]')
xlim([datetime(2015,6,13,2,0,0), datetime(2015,6,13,6,0,0)])


xlabel('time')


% plot(dts2, )
print2(gcf,'~/figures/evap/KOKCdecomp.png', 'quality', '-r300')

%% load radar temporal profiles
ref = ncread('/h/eol/nbarron/workshop/apar-scripts/evapModel/KTLXProfiles.nc', 'refProfile');
zdr = ncread('/h/eol/nbarron/workshop/apar-scripts/evapModel/KTLXProfiles.nc', 'zdrProfile');
% zdr(isnan(zdr)) = 0;
TTobs = ncread('/h/eol/nbarron/workshop/apar-scripts/evapModel/KTLXProfiles.nc', 'time');
TTobs = datetime(1970,1,1) + seconds(TTobs);
ZTobs = ncread('/h/eol/nbarron/workshop/apar-scripts/evapModel/KTLXProfiles.nc', 'height');
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
    print2(fig, strrep(sprintf('~/figures/evap/KOKC_%s.png', titles{i}), ' ', '_'))
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
print2(gcf,'~/figures/evap/KOKCprofiles.png', 'quality', '-r300')

% movmean(ref(8,:),3)

idxStart = 35
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
tab = readmatrix('/h/eol/nbarron/export/soundings/KOUN-2015-06-13_00:00.txt');

height = tab(:,15);
agl = height - height(1);

temp = tab(:,3);
pres = tab(:,2);
dewpt = tab(:,4);
relh = tab(:,5);


fig = figure('units','inches','position',[0,0,3.25,3]);
skewt(pres, temp, dewpt);
print2('~/figures/evap/KOKCsounding.png')

%% set up model run
ev = evap;
ev.xgrid = [50];
ev.ygrid = [50];
ev.zgrid = 50:50:2500;

es = @(T) 2.53e9*exp(-5420./(T)); %hPa

ev.T = interp1(agl, temp+273.15, ev.zgrid);
ev.p = interp1(agl, pres, ev.zgrid);
ev.pv = interp1(agl, es(temp+273.15).*relh/100, ev.zgrid);

ev.u = zeros(size(ev.T));
ev.v = zeros(size(ev.T));
ev.w = zeros(size(ev.T));

ev.T = reshape(ev.T, [numel(ev.xgrid), numel(ev.ygrid), numel(ev.zgrid)]);
ev.p = reshape(ev.p, [numel(ev.xgrid), numel(ev.ygrid), numel(ev.zgrid)]);
ev.pv = reshape(ev.pv, [numel(ev.xgrid), numel(ev.ygrid), numel(ev.zgrid)]);
ev.u = reshape(ev.u, [numel(ev.xgrid), numel(ev.ygrid), numel(ev.zgrid)]);
ev.v = reshape(ev.v, [numel(ev.xgrid), numel(ev.ygrid), numel(ev.zgrid)]);
ev.w = reshape(ev.w, [numel(ev.xgrid), numel(ev.ygrid), numel(ev.zgrid)]);

T0 =ev.T;
dBZi = obsRef(1) + rand(numel(ev.xgrid), numel(ev.ygrid), numel(ev.zgrid))*0;
Zdri = obsZdr(1) + rand(numel(ev.xgrid), numel(ev.ygrid), numel(ev.zgrid))*0;
dBZi(:,:,1:size(dBZi,3)-1) = -inf;
Zdri(:,:,1:size(dBZi,3)-1) = 0;

ev = ev.initFromReflectivity(dBZi, Zdri);
ev.NP = ev.N;

ev.dt = 1;
ev.st = 3.5*3600./ev.dt;
bandName = "S";

sz = numel(ev.zgrid);
%% storm time?
numSteps = ev.st+60/ev.dt;
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


ev = ev.integrate;
dts = [];
dts = TTobs(idxStart);
iObs = 1;
for i = 2:numSteps
    if (i-1)*ev.dt>=obsSeconds(iObs+1) & i*ev.dt<obsSeconds(iObs+2)
        % update the model with the next observation
        fprintf('updating model with observation')
        NTemp0 = ev.N;

        iObs = iObs + 1;
        dbz = zeros(numel(ev.xgrid), numel(ev.ygrid), numel(ev.zgrid));
        dbz(end) = obsRef(iObs);
        zdr = zeros(numel(ev.xgrid), numel(ev.ygrid), numel(ev.zgrid));
        zdr(end) = obsZdr(iObs);
        ev = ev.initFromReflectivity(dbz, zdr);
        ev.N(1,1,1:end-1,:) = NTemp0(1,1,1:end-1,:);
        Ntemp = ev.N;
        Ntemp(:,:,1:end-1,:) = 0;
        ev.NP = Ntemp;


        %% plot the model output
        

        % r = 0.9951 + 0.02510*ev.D - 0.03644*ev.D.^2 + 0.005303*ev.D.^3 - 0.0002492*ev.D.^4;
        fprintf('updated model with observation')
        iObs = iObs + 1;
   
    end

    % if mod(i,600) == 0
    %     % plot the model output
    %     fprintf('plotting model output')
    %     [TT, ZT] = meshgrid((1:numSteps)*ev.dt, ev.zgrid);
    %     variablesToPlot = {T, qv, dBZhh, Zdr, RR, rhohv, kdp};
    %     titles = {'temperature perturbation', 'water vapor saturation', 'Z_{HH}', 'Z_{DR}', 'Rain Rate', 'Correlation Coefficient', 'Specific Differential Phase'};
    %     cLabels = {'K', '', 'dBZ', 'dBZ', 'mm/hr', ' ', 'deg/km'};
    %     clims = {[-8,8], [0.3,1], [-2*10/6,50], [-2,5], [0, 8], [0.99, 1], [0.5e-2, 0.5e-1]};
    %     colormaps = {colortables('newblue'), cbrewer('seq','PuBu',64), colortables('RADAR32'), colortables('hyperspectral'), cbrewer('seq','RdPu',64), flipud(turbo(64)), flipud(autumn(64))+[0,0,0.8]};
    %     plotTimeProfile(TT, ZT,'KOKC', variablesToPlot, titles, cLabels, clims, colormaps)
    % end
    tic
    ev = ev.integrate;
    
    T(i,:) = squeeze(ev.T(1,1,:))-squeeze(T0(1,1,:));
    theta(i,:) = ev.theta(1,1,:);
    p(i,:) = ev.p(1,1,:);
    qv(i,:) = ev.qv(1,1,:);
    pv(i,:) = ev.pv(1,1,:);
    dBZhh(i,:) = ev.Zhh(1,1,:);
    dBZvv(i,:) = ev.Zvv(1,1,:);
    RR(i,:) = ev.RR(1,1,:);
    rhohv(i,:) = ev.rhohv(1,1,:);
    kdp(i,:) = ev.kdp(1,1,:);
    dts(i) = dts(i-1) + seconds(ev.dt);
    toc
    % r = 0.9951 + 0.02510*ev.D - 0.03644*ev.D.^2 + 0.005303*ev.D.^3 - 0.0002492*ev.D.^4;
end
Zdr = dBZhh - dBZvv;



[TT, ZT] = meshgrid((1:numSteps)*ev.dt, ev.zgrid);
variablesToPlot = {T, qv, dBZhh, Zdr, RR, rhohv, kdp};
titles = {'temperature perturbation', 'water vapor saturation', 'Z_{HH}', 'Z_{DR}', 'Rain Rate', 'Correlation Coefficient', 'Specific Differential Phase'};
cLabels = {'K', '', 'dBZ', 'dBZ', 'mm/hr', ' ', 'deg/km'};
clims = {[-12,12], [0.3,1], [-2*10/6,50], [-2,5], [0, 12], [0.99, 1], [0.5e-2, 0.5]};
colormaps = {colortables('newblue'), cbrewer('seq','PuBu',64), colortables('RADAR32'), colortables('hyperspectral'), cbrewer('seq','RdPu',64), flipud(turbo(64)), flipud(autumn(64))+[0,0,0.8]};
plotTimeProfile(dts, ev.zgrid,'KOKC', variablesToPlot, titles, cLabels, clims, colormaps)