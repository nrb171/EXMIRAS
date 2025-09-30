indexToPlot = 3; % 1: Zdr, 2: Zhh_attenuated, 3: VEL
obsPlotBool = 0;
levelOfAOSPRE = 2500;

folder = sprintf('/h/eol/nbarron/workshop/aospre-runs/squall-line-100m/output/%dm/', levelOfAOSPRE);
files = dir(folder);
files = files(~[files.isdir]); % remove directories
files = files(contains({files.name}, '.nc'));
files = sort(string({files.name})); % sort files by name



%% load the RHI scans
inds = 1:100:length(files);
hold on
xs = [];
ys = [];
zs = [];
ds = [];
ts = [];
Zhhs = [];
Zhhms = [];
Zdrs = [];
Vels = [];
lats = [];
lons = [];
files(3:3:18) = []; % remove first 3 files (not RHI scans)
variables = ["Zdr", "Zhh_attenuated", "VEL"];
colormaps = {[winter(12);flipud(autumn(12));flipud(cool(12))], colortables('RADAR32'), colortables('newblue')};
clims = {[-1,5], [1.5*[-2*10/6,50]], [-25,25]}; 
if obsPlotBool
    fig = figure("Units","inches","Position", [0,0,3,3]);
    ax1 = axes(fig);
    hold on
end
for i = 20
    
    file = fullfile(folder, files(i));
    fprintf('Loading %s\n', file)
    
    
    
    x = ncread(file, 'x'); x = x * 100;
    y = ncread(file, 'y'); y = y * 100;
    z = ncread(file, 'z'); 

    d = sqrt((x-x(1,1)).^2 + (y-y(1,1)).^2);
    
    range = ncread(file, 'range');
    dataToPlot = ncread(file, variables(indexToPlot));
    Zdr = ncread(file, 'Zdr');
    Vel = ncread(file, 'VEL');
    Zhh = ncread(file, 'Zhh_attenuated');

    dRegrid = 0:100:max(d(:));
    zRegrid = 0:100:16000;

    ZhhRegrid = 10*log10(griddata(d(:), z(:), 10.^(Zhh(:)/10), dRegrid, zRegrid',"cubic"));
    ZdrRegrid = griddata(d(:), z(:), Zdr(:), dRegrid, zRegrid',"cubic");
    VelRegrid = griddata(d(:), z(:), Vel(:), dRegrid, zRegrid',"cubic");


    figure
    isc = imagesc(dRegrid, zRegrid, real(ZhhRegrid));
    setMeteoColormap(gca, 'Zhh');
    set(gca, 'YDir', 'normal');
    % isc.AlphaData = ones(size(isc.CData));
    % isc.AlphaData(isnan(real(ZhhRegrid'))) = 0;
    print2
    
    figure
    isc = imagesc(dRegrid, zRegrid, real(ZdrRegrid));
    setMeteoColormap(gca, 'Zdr');
    set(gca, 'YDir', 'normal');
    colorbar
    % isc.AlphaData = ones(size(isc.CData));
    % isc.AlphaData(isnan(real(ZdrRegrid'))) = 0;
    print2

    figure
    isc = imagesc(dRegrid, zRegrid, real(VelRegrid));
    setMeteoColormap(gca, 'RadialVelocity');
    set(gca, 'YDir', 'normal');
    % isc.AlphaData = ones(size(isc.CData));
    % isc.AlphaData(isnan(real(VelRegrid'))) = 0;
    print2


    figure
    hold on
    plot(ZhhRegrid(1:21,170), zRegrid(1:21))
    plot(ZdrRegrid(1:21,170), zRegrid(1:21))
    ylim([0,2000])
    print2
end


%{
from netCDF4 import Dataset
from wrf import getvar
import numpy as np

ncfile = Dataset("/export/wind1/bradklotz/CM1_SQUALL_LINE_100m/wrfout_0014420s.nc")

# Get the Sea Level Pressure
rh = getvar(ncfile, "rh")
T = getvar(ncfile, "temp")

np.mean(np.mean(T, axis=1), axis=1)
np.mean(np.mean(P, axis=1), axis=1)

%}
P = [956.0505  , 944.89514 , 934.15216 , 923.49915 , 912.9483  , ...
       902.4854  , 892.13245 , 881.8648  , 871.68286 , 861.60767 , ...
       851.62335 , 841.734   , 831.95    , 822.269   , 812.66174 , ...
       803.14825 , 793.71564 , 784.37335 , 775.11115 , 765.92834 , ...
       756.83264 , 747.81433 , 738.8697  , 730.0115  , 721.2249  , ...
       712.5225  , 703.89484 , 695.34656 , 686.87524 , 678.485   , ...
       670.1784  , 661.9465  , 653.7902  , 645.7213  , 637.7306  , ...
       629.8238  , 621.9907  , 614.24133 , 606.57715 , 598.98083 , ...
       591.46906 , 584.02875 , 576.675   , 569.3864  , 562.18585 , ...
       555.04865 , 547.9888  , 541.0012  , 534.0908  , 527.2473  , ...
       520.47876 , 513.7741  , 507.1445  , 500.5886  , 494.09683 , ...
       487.67343 , 481.3195  , 475.02924 , 468.8109  , 462.6535  , ...
       456.5679  , 450.54083 , 444.58524 , 438.68857 , 432.85666 , ...
       427.089   , 421.3812  , 415.73676 , 410.15366 , 404.6316  , ...
       399.16925 , 393.7678  , 388.42474 , 383.141   , 377.91483 , ...
       372.7467  , 367.63593 , 362.58188 , 357.5844  , 352.64267 , ...
       347.75662 , 342.9238  , 338.14618 , 333.42148 , 328.75043 , ...
       324.13174 , 319.56497 , 315.04855 , 310.58426 , 306.17007 , ...
       301.80624 , 297.49158 , 293.22565 , 289.0082  , 284.83932 , ...
       280.7182  , 276.64288 , 272.6155  , 268.6337  , 264.69824 , ...
       260.80795 , 256.96365 , 253.16339 , 249.40778 , 245.69571 , ...
       242.02766 , 238.40317 , 234.82132 , 231.28297 , 227.78679 , ...
       224.33331 , 220.92058 , 217.54967 , 214.2208  , 210.93224 , ...
       207.68529 , 204.477   , 201.31125 , 198.18391 , 195.09456 , ...
       192.04752 , 189.04083 , 186.07834 , 183.16112 , 180.28806 , ...
       177.45963 , 174.67525 , 171.93365 , 169.23479 , 166.57744 , ...
       163.96094 , 161.385   , 158.84862 , 156.35117 , 153.892   , ...
       151.47075 , 149.08632 , 146.7379  , 144.42564 , 142.14883 , ...
       139.9067  , 137.69803 , 135.52357 , 133.38194 , 131.27293 , ...
       129.19595 , 127.15058 , 125.13655 , 123.15304 , 121.20004 , ...
       119.27665 , 117.383316, 115.51897 , 113.68385 , 111.87726 , ...
       110.098976, 108.348236, 106.62499 , 104.92864 , 103.25827 ];
rh =[70.005554, 71.90741 , 74.13464 , 76.44928 , 78.843094, 81.26029 , ...
       83.696396, 86.16115 , 88.49    , 89.63645 , 87.727325, 84.644485, ...
       82.44525 , 80.813255, 79.324326, 77.92408 , 76.59908 , 75.32339 , ...
       74.08524 , 72.86493 , 71.655365, 70.46351 , 69.28763 , 68.122345, ...
       66.96745 , 65.84524 , 64.77296 , 63.746696, 62.752678, 61.77586 , ...
       60.81679 , 59.87226 , 58.932793, 57.996887, 57.070843, 56.16127 , ...
       55.27646 , 54.423916, 53.593708, 52.790478, 52.019577, 51.275387, ...
       50.553062, 49.851406, 49.17081 , 48.515743, 47.886917, 47.27462 , ...
       46.680313, 46.11076 , 45.56115 , 45.027325, 44.51473 , 44.018436, ...
       43.541283, 43.093273, 42.665367, 42.24522 , 41.845943, 41.478046, ...
       41.134705, 40.81058 , 40.505707, 40.21519 , 39.943676, 39.701027, ...
       39.484653, 39.29155 , 39.13281 , 39.00638 , 38.902657, 38.837776, ...
       38.81838 , 38.8245  , 38.835663, 38.85654 , 38.900204, 38.964672, ...
       39.044712, 39.14085 , 39.256382, 39.39651 , 39.535553, 39.646084, ...
       39.739838, 39.839733, 39.947968, 40.06476 , 40.174026, 40.284283, ...
       40.391163, 40.49499 , 40.59468 , 40.657   , 40.667213, 40.630775, ...
       40.544064, 40.406395, 40.236866, 40.03721 , 39.828094, 39.627254, ...
       39.379807, 39.050343, 38.656723, 38.20484 , 37.68962 , 37.11259 , ...
       36.475624, 35.82982 , 35.197655, 34.54151 , 33.7989  , 32.987354, ...
       32.167892, 31.345547, 30.527958, 29.723314, 28.967255, 28.383598, ...
       28.039577, 27.60882 , 26.84561 , 25.906193, 25.021511, 24.267292, ...
       23.635548, 23.096617, 22.640524, 22.247051, 21.922058, 21.645998, ...
       21.403244, 21.200031, 21.035698, 20.90546 , 20.798489, 20.701181, ...
       20.629314, 20.580713, 20.542852, 20.517916, 20.514359, 20.532207, ...
       20.560818, 20.597746, 20.64423 , 20.693703, 20.737932, 20.767057, ...
       20.768618, 20.72759 , 20.644281, 20.536068, 20.420572, 20.309284, ...
       20.208647, 20.121351, 20.047478, 19.985468];

T = [297.30597, 296.68796, 296.00015, 295.29218, 294.5859 , 293.87408, ...
       293.16592, 292.45074, 291.75427, 291.17065, 290.82468, 290.47232, ...
       289.96375, 289.36215, 288.71317, 288.026  , 287.3035 , 286.55222, ...
       285.77847, 284.98404, 284.17267, 283.35147, 282.52078, 281.68936, ...
       280.85803, 280.03442, 279.2208 , 278.42255, 277.64368, 276.8879 , ...
       276.15973, 275.4618 , 274.7968 , 274.166  , 273.56305, 272.97458, ...
       272.3873 , 271.79745, 271.20447, 270.6079 , 270.00705, 269.40292, ...
       268.79584, 268.18524, 267.57138, 266.95447, 266.33472, 265.71228, ...
       265.088  , 264.46164, 263.83322, 263.203  , 262.56998, 261.93607, ...
       261.2998 , 260.66312, 260.02475, 259.38583, 258.7459 , 258.10507, ...
       257.46362, 256.82217, 256.1805 , 255.53883, 254.89629, 254.25325, ...
       253.60977, 252.96584, 252.32162, 251.67645, 251.03009, 250.38303, ...
       249.73491, 249.0854 , 248.4346 , 247.78171, 247.1265 , 246.46973, ...
       245.80992, 245.14717, 244.48192, 243.8137 , 243.14241, 242.46788, ...
       241.78879, 241.10504, 240.41585, 239.72137, 239.02281, 238.3203 , ...
       237.6145 , 236.90504, 236.1916 , 235.47412, 234.75302, 234.0295 , ...
       233.30379, 232.57602, 231.84546, 231.11421, 230.38374, 229.65408, ...
       228.92616, 228.20088, 227.47804, 226.75896, 226.04483, 225.33463, ...
       224.62796, 223.92528, 223.22621, 222.53209, 221.84375, 221.16093, ...
       220.48518, 219.81668, 219.15381, 218.49542, 217.83873, 217.2015 , ...
       216.66414, 216.32498, 216.16821, 216.09955, 216.05128, 216.00337, ...
       215.95474, 215.90555, 215.85316, 215.79646, 215.73427, 215.66522, ...
       215.58827, 215.50389, 215.41241, 215.31508, 215.21184, 215.10269, ...
       214.98738, 214.86748, 214.74469, 214.61896, 214.49081, 214.36081, ...
       214.22946, 214.09793, 213.96646, 213.83572, 213.70706, 213.58463, ...
       213.4744 , 213.38275, 213.30763, 213.24243, 213.18019, 213.11655, ...
       213.04865, 212.97537, 212.89656, 212.81274];

es = @(T) 2.53e9*exp(-5420./(T)); %hPa
Zm = 50:100:16000

Pv = es(T).*rh/100; % vapor pressure in hPa

temp = griddedInterpolant(Zm, T);
pres = griddedInterpolant(Zm, P); % hPa
pressv = griddedInterpolant(Zm, Pv); % hPa

Rv = 461.5; % J/kg/K
Rd = 287.05; % J/kg/K


%% setup for EXMIRAS runs.
ex = exmiras;
%% set up initial droplet size distribution
ex.xgrid = [50];
ex.ygrid = [50];
ex.zgrid = 25:50:2025;

[ym, xm, zm] = meshgrid(ex.ygrid, ex.xgrid, ex.zgrid);

% es = @(T) 2.53e9*exp(-5420./(T)); %hPa
ex.T = temp(zm);
ex.p = pres(zm);
ex.pv = pressv(zm);

%% setup zhh and zdr time profile
% speed = 15; % m/s
% obsRef = Zhhs(:,scanID);
% obsZdr = Zdrs(:,scanID);
% startInd = find((~isnan(obsRef) | ~isnan(obsZdr) ) & obsZdr > 0.25 & obsRef>10,1);
% startInd = 1;
% obsRef = obsRef(startInd:720);
% obsZdr = obsZdr(startInd:720);
    % timeToGate = (0:75:75*numel(obsRef)-1)/speed;




    ex.u = zeros(size(ex.T));
    ex.v = zeros(size(ex.T));
    ex.w = zeros(size(ex.T));

    %% reshape to match model grid
    ex.T = reshape(ex.T, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);
    ex.p = reshape(ex.p, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);
    ex.pv = reshape(ex.pv, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);
    ex.u = reshape(ex.u, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);
    ex.v = reshape(ex.v, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);
    ex.w = reshape(ex.w, [numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid)]);

    ex.WaterVaporInteraction = false;


    %% initialize dsd from reflectivity/zdr
    dBZi = ZhhRegrid(21, 170) + rand(numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid))*0;
    Zdri = ZdrRegrid(21, 170) + rand(numel(ex.xgrid), numel(ex.ygrid), numel(ex.zgrid))*0;
    dBZi(:,:,1:size(dBZi,3)-1) = -inf;
    Zdri(:,:,1:size(dBZi,3)-1) = 0;


    % initialize dsd from reflectivity/zdr
    ex = ex.initFromLambdaName("C");
    ex = ex.initFromReflectivity(dBZi, Zdri);
    ex.NP = ex.N;

    %set time step/total integration time
    ex.dt = 0.5;
    ex.st = 600/ex.dt;

    sz = numel(ex.zgrid);
    %% storm time?
    numSteps = ex.st;
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
    lvl = 0;
    % figure
    % hold on
    % % plot(ex.D,squeeze(ex.N(1,1,end-lvl,:)), "DisplayName", "N")
    % plot(ex.D, squeeze(ex.dNevap(1,1,end-lvl,:)), "DisplayName", "dNevap")
    % plot(ex.D, squeeze(ex.dNfallout(1,1,end-lvl,:)), "DisplayName", "dNfallout")
    % plot(ex.D, squeeze(ex.dNcoal(1,1,end-lvl,:)),"DisplayName", "dNcoal")
    % % plot(ex.D, squeeze(ex.N(1,1,end-lvl,:)) + squeeze(ex.dNevap(1,1,end-lvl,:)) + squeeze(ex.dNfallout(1,1,end-lvl,:)) + squeeze(ex.dNcoal(1,1,end-lvl,:)), "DisplayName", "N(t+dt)")
    % % plot(obj.D, squeeze(dDevap(1,1,end-1,:)), "DisplayName", "dN")
    % legend
    % print2(gcf, '.temp.png')

    % dts = [];
    % dts = TTobs(1);
    iObs = 1;
    % fig = figure('units','inches','position',[0,0,6.05,3]);
    % hold on
    errorLog = [];
    for i = 2:numSteps
        % try
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
        % dts(i) = dts(i-1) + seconds(ex.dt);

        error = rms(interp1(ex.zgrid, squeeze(ex.Zhh(1,1,:)), zRegrid(1:21)) - ZhhRegrid(1:21,170)', 'omitnan');
        error2 = rms(interp1(ex.zgrid, squeeze(ex.Zdr(1,1,:)), zRegrid(1:21)) - ZdrRegrid(1:21,170)', 'omitnan');  
        errorLog(end+1) = mean([error, error2]);
        fprintf('Step %d/%d, RMS error: %.2f\n', i, numSteps, errorLog(end));

        if i > 100
            if mod(i,50) == 0
                figure
                hold on
                plot(squeeze(ex.Zhh(1,1,:)), ex.zgrid, 'b')
                plot(squeeze(ex.Zdr(1,1,:)), ex.zgrid, 'r')
                plot(ZhhRegrid(1:21,170), zRegrid(1:21), 'c--')
                plot(ZdrRegrid(1:21,170), zRegrid(1:21), 'm--')
                legend('Zhh sim', 'Zdr sim', 'Zhh obs', 'Zdr obs')
                print2
            end
            % errorLog(end+1) = mean([error, error2]);

            errorTendency = movmean(gradient(errorLog), 10);

            if all(errorTendency(end-10:end)>0)
                fprintf('Error increasing, stopping early at step %d/%d\n', i, numSteps)
                figure
                hold on
                plot(squeeze(ex.Zhh(1,1,:)), ex.zgrid, 'b')
                plot(squeeze(ex.Zdr(1,1,:)), ex.zgrid, 'r')
                plot(ZhhRegrid(1:21,170), zRegrid(1:21), 'c--')
                plot(ZdrRegrid(1:21,170), zRegrid(1:21), 'm--')
                legend('Zhh sim', 'Zdr sim', 'Zhh obs', 'Zdr obs')
                print2
                break
            end


        end
    end
    Zdr = dBZhh - dBZvv;

    save(sprintf('./SL-%d_%d.mat', levelOfAOSPRE,scanID), 'T', 'p', 'qv', 'pv', 'dBZhh', 'dBZvv', 'Zdr', 'RR', 'rhohv', 'kdp', 'ex');
