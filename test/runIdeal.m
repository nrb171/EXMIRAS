%%! set up initial droplet size distribution
ev = evap;
%% set up initial droplet size distribution
ev.xgrid = [50];
ev.ygrid = [50];
ev.zgrid = 50:50:2050;
bandName = "S"

sx = numel(ev.xgrid);
sy = numel(ev.ygrid);
sz = numel(ev.zgrid);
zdr = 0.5
dBZi = 30 + rand(sx, sy, sz)*0;
Zdri = zdr + rand(sx, sy, sz)*0;
dBZi(:,:,1:size(dBZi,3)-1) = -inf;
Zdri(:,:,1:size(dBZi,3)-1) = 0;

ev = ev.initFromLambdaName("S");
ev = ev.initFromReflectivity(dBZi, Zdri);

ev.NP = ev.N;

%%! set up initial state variables
temp = @(z) -6.5*z/1000 + 293;
pres = @(z) 1000*exp(-z/1000/8.4);
es = @(T) 2.53e9*exp(-5420./(T)); %hPa
[ym, xm, zm] = meshgrid(ev.ygrid, ev.xgrid, ev.zgrid);
ev.T = temp(zm);
ev.p = pres(zm);
ev.pv = es(ev.T)*0.5;

ev.u = zeros(size(ev.N, [1:3]))
ev.v = zeros(size(ev.N, [1:3]))
ev.w = zeros(size(ev.N, [1:3]))


plot(ev.D, squeeze(ev.N(1,1,end,:)) )
print2(gcf, './.temp.png')


ev.dt = 0.5;
ev.st = 60./ev.dt; 

fig = figure(units="inches", position=[0,0,3.75,3]);
plot(ev.D, (squeeze(ev.dNcoal(1,1,end,:))) )
hold on
plot(ev.D, (squeeze(ev.dNevap(1,1,end,:)) ))
yyaxis right
plot(ev.D, squeeze(ev.dNfallout(1,1,end,:)))

% plot(ev.D, (squeeze(ev.dNcoal(1,1,end,:)) )+squeeze(ev.dNevap(1,1,end,:)) )
yyaxis left
xscale("log")
xlabel('D [mm]')
ylabel('dN [m^{-3} mm^{-1}]')
legend('coal', 'evap', 'fallout')
print2(fig, './.temp.png')



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
profile on
for i = 1:numSteps
    tic
    ev = ev.integrate;
    T(i,:) = squeeze(ev.T(1,1,:))-squeeze(temp(zm(1,1,:)));
    theta(i,:) = ev.theta(1,1,:);
    p(i,:) = ev.p(1,1,:);
    qv(i,:) = ev.qv(1,1,:);
    pv(i,:) = ev.pv(1,1,:);
    dBZhh(i,:) = ev.Zhh(1,1,:);
    dBZvv(i,:) = ev.Zvv(1,1,:);
    RR(i,:) = ev.RR(1,1,:);
    rhohv(i,:) = ev.rhohv(1,1,:);
    kdp(i,:) = ev.kdp(1,1,:);
    toc
    % r = 0.9951 + 0.02510*ev.D - 0.03644*ev.D.^2 + 0.005303*ev.D.^3 - 0.0002492*ev.D.^4;
end
profsave
Zdr = dBZhh - dBZvv;




[TT, ZT] = meshgrid((1:numSteps)*ev.dt, ev.zgrid);



variablesToPlot = {T, qv, dBZhh, Zdr, RR, rhohv, kdp};
titles = {'temperature perturbation', 'water vapor saturation', 'Z_{HH}', 'Z_{DR}', 'Rain Rate', 'Correlation Coefficient', 'Specific Differential Phase'};
cLabels = {'K', '', 'dBZ', 'dBZ', 'mm/hr', ' ', 'deg/km'};
clims = {[-8,8], [0.3,1], [-2*10/6,50], [-2,5], [0, 8], [0.99, 1], [0.5e-2, 0.5e-1]};
colormaps = {colortables('newblue'), cbrewer('seq','PuBu',64), colortables('RADAR32'), colortables('hyperspectral'), cbrewer('seq','RdPu',64), flipud(turbo(64)), flipud(autumn(64))+[0,0,0.8]};
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
    print2(fig, strrep(sprintf('~/figures/evap/%s_ideal_%s.png', bandName,titles{i}), ' ', '_'))
end