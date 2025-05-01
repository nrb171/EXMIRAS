
%% call exmiras
ex = exmiras;
%% USER: set simulation grid resolution
ex.xgrid = %!xgrid in meters, e.g. [50] Ideally, don't use more than 1 for horizontal grids
ex.ygrid = %!ygrid in meters, e.g. [50]
ex.zgrid = %!zgrid in meters, e.g. 50:50:2050. These should be set so that 0.5*ex.dt.^2*9.81 - ex.w*ex.dt << ex.zgrid(2)-ex.zgrid(1) to avoid numerical instabilities. 

%% USER: set initial conditions
bandName = %!Letter of the band, e.g. "S"
dBZStart = %!dBZ at the top of the simulation, e.g. 20
zdr = %!Zdr at the top of the simulation, e.g. 0.5
T0m = %!starting Temperature at the bottom of the simulation, e.g. 273.15
lapseRate = %!lapse rate in K/m, e.g. -6.5
humidity = %!relative humidity in decimal, e.g. 0.3

%% USER: set simulation time
ex.dt = %!time step in seconds, e.g. 0.5;
ex.st = %!simulation time in number of steps, e.g. 3*3600./ex.dt;
humidityLimit = %!limit for relative humidity in decimal, e.g. 0.95. If all values in the prior step of the simulation are 1 it will stop the simulation.

%% USER: assign variables to be saved
variablesToSave = %!variables to save, e.g. {'T', 'p', 'pv', 'qv', 'Zhh', 'Zvv', 'RR', 'rhohv', 'kdp', 'Zdr', 'theta'};


%% Additional settings: 
% - You can add vertical motion to the simulation by setting ex.w to something other than zero.

%%%%%%%%%%%%%%%%%%%%%% END USER INPUT %%%%%%%%%%%%%%%%%%%%%%

%% set up the initial droplet size distribution
sx = numel(ex.xgrid);
sy = numel(ex.ygrid);
sz = numel(ex.zgrid);
dBZi = dBZStart + rand(sx, sy, sz)*0;
Zdri = zdr + rand(sx, sy, sz)*0;
dBZi(:,:,1:size(dBZi,3)-1) = -inf;
Zdri(:,:,1:size(dBZi,3)-1) = 0;

% run the dsd solver
ex = ex.initFromLambdaName(bandName);
ex = ex.initFromReflectivity(dBZi, Zdri);
ex.NP = ex.N;

%% set up initial state variables
temp = @(z) lapseRate*z/1000 + T0m;
pres = @(z) 1000*exp(-z/1000/8.4);
es = @(T) 2.53e9*exp(-5420./(T)); %hPa
[ym, xm, zm] = meshgrid(ex.ygrid, ex.xgrid, ex.zgrid);
ex.T = temp(zm);
ex.p = pres(zm);
ex.pv = es(ex.T)*humidity;

%% set up initial kinematic variables
ex.u = zeros(size(ex.N, [1:3]))
ex.v = zeros(size(ex.N, [1:3]))
ex.w = zeros(size(ex.N, [1:3]))


%% set up structure to save the simulation
ExmirasRun = struct();
ExmirasRun.InitVariables.p = pres(zm);
ExmirasRun.InitVariables.T = temp(zm);
ExmirasRun.InitVariables.pv = es(ex.T)*humidity;
ExmirasRun.ID = sprintf('%s_ideal_%d_%1.2f_%2.0d', bandName, T0m, humidity,dBZStart);

for i = 1:numel(variablesToSave)
    %note: variablesToSave is set above by the user
    ExmirasRun.(variablesToSave{i}) = NaN(ex.st, sx, sy, sz);
end

%% run the simulation
for i = 1:numSteps
    % update user on progress
    if mod(i/ex.st, 0.05) == 0
        fprintf('%.2f%%\n', i/numSteps*100)
    end
    ex = ex.integrate;

    % save variables
    for j = 1:numel(variablesToSave)
        ExmirasRun.(variablesToSave{j})(i,:) = ex.(variablesToSave{j})(1,1,:);
    end
    ExmirasRun.Zdr(i,:) = 10^(ex.Zhh./10) - 10^(ex.Zvv./10);

    % stop integrating if the water vapor saturation is reached above 99%
    if all(ex.qv(1,1,:)>=humidityLimit) | all(ex.qv(1,1,:)>=1)
        break
    end 
end

%% wrap up and save the simulation
ExmirasRun.Grids.tgrid = ((0:i-1))*ex.dt;
ExmirasRun.Grids.zgrid = ex.zgrid;
save(['./',ExmirasRun.ID,'.mat'], 'ExmirasRun')
