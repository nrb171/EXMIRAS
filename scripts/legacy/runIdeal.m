% if you acquired this file from the APAR folder, please see https://github.com/nrb171/EXMIRAS for the latest version of this and other scripts.

% use the job scheduler (/glade/u/home/nbarron/workshop/EXMIRAS/scripts/helpers/submitEXMIRASJobs.sh) to run the ideal simulation on Derecho. This should be a fairly good example of how to run this on other HPC systems.
% T0 = str2num(getenv('T0'));
% humidity = str2num(getenv('HUMIDITY'));
global saveDir
saveDir = '/h/eol/nbarron/work/exmirasData';

tRange = 300;
tTop = numel(tRange);

% runEXMIRAS(300, 0.3, 50, 2, "C", -6.5, [saveDir, 'DSD-Test-fallout'], 'coalToggle', false, 'evapToggle', false, 'waterVaporInteraction', false)
% runEXMIRAS(300, 0.3, 50, 2, "C", -6.5, [saveDir, 'DSD-Test-fallout+evap'], 'coalToggle', false, 'waterVaporInteraction', false)
% runEXMIRAS(300, 0.3, 50, 2, "C", -6.5, [saveDir, 'DSD-Test-fallout+coal'],'evapToggle', false, 'waterVaporInteraction', false)
% runEXMIRAS(300, 0.3, 50, 2, "C", -6.5, [saveDir, 'DSD-Test-all'], 'waterVaporInteraction', false)

runEXMIRAS(300, 0.3, 50, 1.5, "C", -6.5, [saveDir, 'diff-kernel'], 'waterVaporInteraction', false)
runEXMIRAS(300, 0.3, 50, 1, "C", -6.5, [saveDir, 'diff-kernel'], 'waterVaporInteraction', false)
runEXMIRAS(300, 0.3, 50, 0.5, "C", -6.5, [saveDir, 'diff-kernel'], 'waterVaporInteraction', false)
runEXMIRAS(300, 0.3, 50, 2, "C", -6.5, [saveDir, 'diff-kernel'], 'waterVaporInteraction', false)


return

for ii = 1:tTop
    T0 = tRange(ii);
    for humidity = 0.3:0.3:0.9
        for dBZStart = 30:5:50
            for zdr = 0.5:0.5:1.5
                runEXMIRAS(T0, humidity, dBZStart, zdr, "S", -6.5, saveDir)
            end
        end
    end
end


files = dir('/h/eol/nbarron/work');
files = {files.name};
files = files(contains(files, '_ideal_'));

variableToPlot = {'dTdt', 'RH', 'Zhh', 'Zdr'};


for jj = 1:numel(variableToPlot)
    legendStr="";
  
    var2=NaN(numel(files), 7200);
    for ii = 1:numel(files)
        f = load(['/h/eol/nbarron/work','/',files{ii}]);

        switch variableToPlot{jj}
            case 'dTdt'
                [~,dTdtTemp]=gradient(f.ExmirasRun.T, 0.5);
                dTdt = NaN([7200,41]);
                dTdt(1:size(dTdtTemp,1),:) = dTdtTemp*3600; % convert to K/hr
                
                f.ExmirasRun.dTdt = dlog10(dTdt);%sign(dTdt).*sqrt(abs(dTdt));
                

            case 'RH'
                f.ExmirasRun.RH = f.ExmirasRun.qv*100;
        end
        
        legendStr(end+1) = string(files{ii}(12:end-4));
        
        var = f.ExmirasRun.(variableToPlot{jj});
        var2(ii,1:size(var,1)) = var(:,1);

    end


    fig = figure('Units', 'inches', 'Position', [0,0,2.85,3]);
    hold on
    % si=[1,7,13,19,25,31,37,43,49];
    % ei=si+5
    li = 1;
    yticksSave = [];
    for ii = 1:15
        plottedInds = [1*(ii-1)+3*(ii-1):1*(ii-1)+3*(ii-1)+2]+1;
        yticksSave = [yticksSave, plottedInds]; 
        isc=imagesc(linspace(0,3600, 7200),plottedInds, var2(li:li+2,:));

        yline([plottedInds(2:end)-1]+0.5, 'LineWidth', 0.1)
        isc.AlphaData = ones(size(isc.CData));
        isc.AlphaData(isnan(var2(li:li+2,:))) = 0;
        li = li+3;

    end

    setMeteoColormap(gca, variableToPlot{jj})
    set(gca, 'YDir', 'normal')
    % hold on
    xline([600:600:3000], 'LineWidth', 0.5,'Alpha', 0.2)

    % set(gca, 'ColorScale', 'linear')
    xticks(0:600:3600)
    xticklabels(0:10:60)
    yticks(yticksSave)
    yticklabels([])
    xlim([0,3600])
    ylim([0.5, max(yticksSave)+0.5])
    shading flat
    print2(fig, sprintf('~/figures/publications/exmiras/ideal_all_%s.pdf', variableToPlot{jj}), 'quality', '-r300')
    writelines(legendStr, sprintf('~/figures/publications/exmiras/ideal_%s_legend.txt', variableToPlot{jj}))
end


files = dir('/h/eol/nbarron/work');
files = {files.name};
files = files(contains(files, 'DSD-Test') & contains(files, 'C_ideal') & contains(files, "50"));
% files = flip(files);

De = logspace(log10(0.1),log10(8), 250+1);
D =  (De(1:end-1) + De(2:end))/2;

N = [];

ex = exmiras
times = [3, 6, 9];
for kk = 1:3
    yticksSave = [];
    yticklabelsSave = string();
    figure('Units', 'inches', 'Position', [0,0,2,3]);
    hold on
    for jj = 1:numel(files)
        f = load(['/h/eol/nbarron/work','/',files{jj}]);
        N(1,:) = interp1(D, squeeze(f.ExmirasRun.N(times(kk),1,:)), linspace(min(D), max(D), 250));
        N(2,:) = interp1(D, squeeze(f.ExmirasRun.N(times(kk),11,:)), linspace(min(D), max(D), 250));
        N(3,:) = interp1(D, squeeze(f.ExmirasRun.N(times(kk),21,:)), linspace(min(D), max(D), 250));
        N(4,:) = interp1(D, squeeze(f.ExmirasRun.N(times(kk),31,:)), linspace(min(D), max(D), 250));
        N(5,:) = interp1(D, squeeze(f.ExmirasRun.N(times(kk),41,:)), linspace(min(D), max(D), 250));

        N=N./max(N,[],2)*100;


        ex = ex.initFromLambdaName("C");
        ex = ex.initFromReflectivity(40, 1);
        Zhh = [];
        Zdr = [];
        for ii = 1:5
            ex.N(1,1,1,:) = squeeze(f.ExmirasRun.N(times(kk),1+(ii-1)*10,:));
            Zhh(end+1) = ex.Zhh(1,1,1);
            Zdr(end+1) = ex.Zdr(1,1,1);
        end
        ex.N(:,:,1,:)=squeeze(f.ExmirasRun.N(times(kk),1,:));
        % Zhh(1) = f.ExmirasRun.Zhh(times(kk)*120,1);
        % Zhh(2) = f.ExmirasRun.Zhh(times(kk)*120,11);
        % Zhh(3) = f.ExmirasRun.Zhh(times(kk)*120,21);
        % Zhh(4) = f.ExmirasRun.Zhh(times(kk)*120,31);
        % Zhh(5) = f.ExmirasRun.Zhh(times(kk)*120,41);


        N2 = squeeze(f.ExmirasRun.N(times(kk),1:10:41,:));
        % N2(2,:) = interp1(D, squeeze(f.ExmirasRun.N(times(kk),11,:)), D);
        % N2(3,:) = interp1(D, squeeze(f.ExmirasRun.N(times(kk),21,:)), D);
        % N2(4,:) = interp1(D, squeeze(f.ExmirasRun.N(times(kk),31,:)), D);
        % N2(5,:) = interp1(D, squeeze(f.ExmirasRun.N(times(kk),41,:)), D);

        

        log10(1+N);
        if jj>1
            % N = N-N2;
            imagesc(D, [1:5]+5*(jj-1)+1*(jj-1), N, 'AlphaData', ~isnan(log10(1+N)))
            yticksSave = [yticksSave, [1:5]+5*(jj-1) + 1*(jj-1)-0.5];

            [~, maxInds] = max(N2, [], 2);
            means = sum(N2.*D,2)./sum(N2,2);

            for ii = 1:numel(maxInds)
                yinds = [1:5]+5*(jj-1) + 1*(jj-1);
                plot(D([maxInds(ii), maxInds(ii)]), [yinds(ii)-0.5, yinds(ii)+0.5], 'w', 'LineWidth', 0.5)
                plot([means(ii), means(ii)], [yinds(ii)-0.5, yinds(ii)+0.5], 'r', 'LineWidth', 0.5)
                text(3.25, yinds(ii),string(round(Zhh(ii),2)), "FontSize", 3)
                text(4.1, yinds(ii),string(round(Zdr(ii),2)), "FontSize", 3)
                
            end

            
        else
            % N2(1,:) = squeeze(f.ExmirasRun.N(1,41,:));
            % N2(2,:) = squeeze(f.ExmirasRun.N(1,41,:));
            % N2(3,:) = squeeze(f.ExmirasRun.N(1,41,:));
            % N2(4,:) = squeeze(f.ExmirasRun.N(1,41,:));
            % N2(5,:) = squeeze(f.ExmirasRun.N(1,41,:));
            imagesc(D, [1:5]+5*(jj-1), (N), 'AlphaData', ~isnan(log10(1+N)));
            yticksSave = [yticksSave, [1:5]+5*(jj-1)-0.5];

            [~, maxInds] = max(N2, [], 2);
            means = sum(N2.*D,2)./sum(N2,2);
            for ii = 1:numel(maxInds)
                yinds = [1:5]+5*(jj-1);
                plot(D([maxInds(ii), maxInds(ii)]), [yinds(ii)-0.5, yinds(ii)+0.5], 'w', 'LineWidth', 0.5)
                plot([means(ii), means(ii)], [yinds(ii)-0.5, yinds(ii)+0.5], 'r', 'LineWidth', 0.5)
                text(3.25, yinds(ii), string(round(Zhh(ii),2)), "FontSize", 3)
                text(4.1, yinds(ii), string(round(Zdr(ii),2)), "FontSize", 3)
                
            end

            

        end
        set(gca, 'YDir', 'normal')
        yticklabelsSave = [yticklabelsSave, string([1:10:41]*50 - 50)];
    end
    clim([0,100])
    yticks(yticksSave+0.5)
    ylim([0, max(yticksSave)+1])
    grid on
    xticks([0:0.5:8])
    xtlabels = string([0:0.5:8]);
    xtlabels(2:2:end) = "";
    xticklabels(xtlabels)
    xtickangle(0)
    % set(gca, 'XScale','log')
    xlim([0,4.5])
    yticklabels(yticklabelsSave(2:end))
    colormap(cbrewer('seq', 'YlGnBu', 64))
    print2(gcf, sprintf('~/figures/publications/exmiras/ideal_DSD_vertical_%d.pdf', kk), 'quality', '-r300')

    %% generate colorbar
    fig = figure('Units', 'inches', 'Position', [0,0,0.5,2]);
    hold on
    colormap(cbrewer('seq', 'YlGnBu', 64))
    caxis([0,150])
    colorbar
    ax = gca;
    ax.Visible = 'off';
    print2(fig, sprintf('~/figures/publications/exmiras/ideal_DSD_vertical_colorbar.pdf'), 'quality', '-r300')
    
end

f = load(['/h/eol/nbarron/work','/',files{end}]);

%% plot dN at the start of the simulation
fig = figure('Units', 'inches', 'Position', [0,0,2,2]);
hold on
plot(D, (squeeze(f.ExmirasRun.dNevap(1,end,:))))
plot(D, (squeeze(f.ExmirasRun.dNcoal(1,end,:))))
plot(D, (squeeze(f.ExmirasRun.dNcoal(1,end,:))) + (squeeze(f.ExmirasRun.dNevap(3,end,:))), '--k')
ylim(([-0.5, 0.5]))
xlim([0,4])
grid on
% yticks(([-0.05:0.025:0.05]))
% yticklabels([-10:2:10])
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
yyaxis right
plot(D, (squeeze(f.ExmirasRun.dNfallout(3,end,:))))
ylim([-50, 50])
xlim([0,4])
% plot(D, squeeze(f.ExmirasRun.N(10,40,:)))
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
print2(gcf, sprintf('~/figures/publications/exmiras/ideal_DSD_dNs.pdf'))

%% plot dn at the end of the simulation, after steady state
fig = figure('Units', 'inches', 'Position', [0,0,2,2]);
hold on
plot(D, (squeeze(f.ExmirasRun.dNevap(end,21,:))))
plot(D, (squeeze(f.ExmirasRun.dNcoal(end,21,:))))
plot(D, (squeeze(f.ExmirasRun.dNcoal(end,21,:))) + (squeeze(f.ExmirasRun.dNevap(end,21,:))), '--k')
grid on

ylim(([-0.5, 0.5]))
xlim([0,4])
% yticks(([-0.05:0.025:0.05]))
% yticklabels([-10:2:10])
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
yyaxis right
plot(D, (squeeze(f.ExmirasRun.dNfallout(end,21,:))))
ylim([-0.5, 0.5])
xlim([0,4])
% plot(D, squeeze(f.ExmirasRun.N(10,40,:)))
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
print2(gcf, sprintf('~/figures/publications/exmiras/ideal_DSD_dNe.pdf'))

%% get some DSDs to plot up
ex = exmiras;

gamma = [];
mu = [];
zdr1 =[];
N0 = [];
N = [];
colors = winter(100);
zdrs = 1:0.5:3;
zhh=30;
for ii = 1:numel(zdrs)
    zdr = zdrs(ii);
    ex = ex.initFromLambdaName("C");
    ex = ex.initFromReflectivity(zhh, zdr);


    gamma(end+1) = ex.gamma;
    mu(end+1) = ex.mu;
    zdr1(end+1) = zdr;
    N0(end+1) = ex.N00;
    N(ii,:) = ex.N(1,1,end,:);
end
% gammaDist = @
fig = figure("Units","inches","Position", [0,0,2,2]);
hold on;
plot(ex.D, N(1,:));
plot(ex.D, N(2,:));
plot(ex.D, N(3,:));
plot(ex.D, N(4,:));

% plot(ex.D, N(5,:));
set(gca, 'YScale', 'log')
print2(gcf, '~/figures/publications/exmiras/ideal_DSD_diffzdrs.pdf', 'quality', '-r300')



%% calculate relative dN


files = dir('/h/eol/nbarron/work');
files = {files.name};
files = files(contains(files, 'DSD-Test') & contains(files, 'C_ideal') & contains(files, "50"));
% files = 
files = flip(files);
ffallout = load(['/h/eol/nbarron/work','/',files{1}]);
fevap = load(['/h/eol/nbarron/work','/',files{2}]);
fcoal = load(['/h/eol/nbarron/work','/',files{3}]);
f = load(['/h/eol/nbarron/work','/',files{4}]);

Ntop = squeeze(f.ExmirasRun.N(20,end,:));
Nbot = squeeze(f.ExmirasRun.N(20,1,:));
Ntopfallout = squeeze(ffallout.ExmirasRun.N(20,end,:));
Nbotfallout = squeeze(ffallout.ExmirasRun.N(20,1,:));
Ntopcoal = squeeze(fcoal.ExmirasRun.N(20,end,:));
Nbotcoal = squeeze(fcoal.ExmirasRun.N(20,1,:));
Ntopevap = squeeze(fevap.ExmirasRun.N(20,end,:));
Nbotevap = squeeze(fevap.ExmirasRun.N(20,1,:));

dnzfallout=(squeeze(ffallout.ExmirasRun.N(20,1,:))-squeeze(ffallout.ExmirasRun.N(20,end,:)));
dnzcoal=(squeeze(fcoal.ExmirasRun.N(20,1,:))-squeeze(fcoal.ExmirasRun.N(20,end,:)));
dnzevap=(squeeze(fevap.ExmirasRun.N(20,1,:))-squeeze(fevap.ExmirasRun.N(20,end,:)));
De = logspace(log10(0.1),log10(8), 250+1);
D =  (De(1:end-1) + De(2:end))/2;

recon =squeeze(f.ExmirasRun.N(20,end,:))+dnzevap+ dnzcoal;
% normalizer = trapz(D, Nbot)./trapz(D, recon);
% recon = recon*normalizer;

figure('Units', 'inches', 'Position', [0,0,3,3]);
% plot(D,dnzevap, 'DisplayName', 'evap')
hold on
% plot(D,dnzcoal, 'DisplayName', 'coal')
plot(D, Nbot-Ntop, 'DisplayName', 'diff')

% func = @(N,mu,lambda) rms(-N.*D.^mu.*exp(-lambda.*D)-(Nbot-Ntop)');
gamma = @(N,mu,lambda) N.*D.^mu.*exp(-lambda.*D);
funcTop = @(N,mu,lambda) rms(gamma(N,mu,lambda)-(Ntop)');
foptTop = fminsearchbnd(@(x) func(x(1), x(2), x(3)), [2.7408e+04,4, 3.5], [1e4, -2, 0], [4e4,15, 20])
funcBot = @(N,mu,lambda) rms(gamma(N,mu,lambda)-(Nbot)');
foptBot = fminsearchbnd(@(x) funcBot(x(1), x(2), x(3)), [2.7408e+04,4, 3.5], [1e3, -2, 0], [4e4,15, 20])
funcDiff = @(N,mu,lambda) rms(gamma(N,mu,lambda)-((Nbot-Ntop)'));
fopt = fminsearchbnd(@(x) funcDiff(-x(1), x(2), x(3)), foptTop, [foptTop(1)-1, foptTop(2)-0.01, 0], [foptTop(1)+1, foptTop(2)+0.01, 20])
% figure

foptTop(1), foptBot(1), fopt(1)
foptTop(2), foptBot(2), fopt(2)
foptTop(3)-fopt(3)

figure
% plot(D, gamma(2.3e+04,4, 3.5))
hold on
% plot(D, gamma(-fopt(1),fopt(2), fopt(3)), '--')
% plot(D, Nbot-Ntop)
plot(D, Nbot)
plot(D, Ntop-gamma(fopt(1), fopt(2), fopt(3)), '--')

yyaxis right
plot(D, gamma(fopt(1), fopt(2), fopt(3))./gamma(foptTop(1),foptTop(2), foptTop(3)))
hold on
plot(D,exp(-(fopt(3)-foptTop(3))*D))
% plot(D, Ntop)
print2
% plot(D, -fopt(1)*D.^fopt(2).*exp(-fopt(3)*D), '--', 'DisplayName', 'diff (fit guess)')
% print2
% pf = polyfit(D, Nbot-Ntop, 5);
% plot(D, polyval(pf, D), '--', 'DisplayName', 'diff (fit)')
% plot(D,dnzfallout, 'DisplayName', 'fallout')

% yyaxis right
% plot(D, squeeze(f.ExmirasRun.N(20,1,:)),"DisplayName", "bottom")
hold on
plot(D, Nbot, "DisplayName", "bottom")
% plot(D, Ntopcoal+dnzcoal*0.6+dnzevap*0.6, "DisplayName", "bottom (calc)")
plot(D, Ntop, "DisplayName", "top")
% yscale('log')
% plot(D, squeeze(f.ExmirasRun.N(20,end,:)) + dnzevap'+dnzcoal', "DisplayName", "diff")
legend


yyaxis right
plot(D,(-fopt(1)*D.^fopt(2).*exp(-fopt(3)*D))'./Ntop)
print2(gcf, '~/figures/publications/exmiras/ideal_DSD_dN_diff.png', 'quality', '-r300')



%%% function to run the ideal simulation %%%
function runEXMIRAS(T0m, humidity, dBZStart, zdr, bandName, lapseRate, saveDir, varargin)
    tic
    % Run the ideal simulation                  
    %%! set up initial droplet size distribution
    fprintf('T0m = %1.2f, humidity = %1.2f, dBZStart = %d, zdr = %1.2f\n', T0m, humidity, dBZStart, zdr)
    ex = exmiras;
    %% set up initial droplet size distribution
    ex.xgrid = [50];
    ex.ygrid = [50];
    ex.zgrid = 50:50:2050;

    sx = numel(ex.xgrid);
    sy = numel(ex.ygrid);
    sz = numel(ex.zgrid);
    dBZi = dBZStart + rand(sx, sy, sz)*0;
    Zdri = zdr + rand(sx, sy, sz)*0;
    dBZi(:,:,1:size(dBZi,3)-1) = -inf;
    Zdri(:,:,1:size(dBZi,3)-1) = 0;

    ex = ex.initFromLambdaName(bandName);
    ex = ex.initFromReflectivity(dBZi, Zdri);

    ex.NP = ex.N;

    %%! set up initial state variables
    temp = @(z) lapseRate*z/1000 + T0m;
    pres = @(z) 1000*exp(-z/1000/8.4);
    es = @(T) 2.53e9*exp(-5420./(T)); %hPa
    [ym, xm, zm] = meshgrid(ex.ygrid, ex.xgrid, ex.zgrid);
    ex.T = temp(zm);
    ex.p = pres(zm);
    ex.pv = es(ex.T)*humidity;

    ex.u = zeros(size(ex.N, [1:3]));
    ex.v = zeros(size(ex.N, [1:3]));
    ex.w = zeros(size(ex.N, [1:3]));

    ex.dt = 0.5;
    % ex.st = 3600./ex.dt; 
    ex.st = 1200./ex.dt; 
    


    %% set storm time
    numSteps = ex.st+0;
    variablesToSave = {'T', 'p', 'pv', 'qv', 'Zhh', 'Zvv', 'RR', 'rhohv', 'kdp', 'Zdr', 'theta'};
    ExmirasRun = struct();
    ExmirasRun.initVariables.p = pres(zm);
    ExmirasRun.initVariables.T = temp(zm);
    ExmirasRun.initVariables.pv = es(ex.T)*humidity;
    ExmirasRun.ID = sprintf('%s_ideal_%d_%1.2f_%2.0d_%1.2f', bandName, T0m, humidity,dBZStart, zdr);

    %% check to see if the file already exists, if so, break the loop
    filesInSaveDir = dir(saveDir);
    if any(contains({filesInSaveDir.name}, ExmirasRun.ID))
        fprintf('File %s already exists, skipping...\n', ExmirasRun.ID)
        return
    end

    %% set up the struct to save the run
    for i = 1:numel(variablesToSave)
        ExmirasRun.(variablesToSave{i}) = NaN(numSteps, sz);
    end
    ExmirasRun.N = [];
    ExmirasRun.dNevap = [];
    ExmirasRun.dNcoal = [];
    ExmirasRun.dNfallout = [];

    if numel(varargin)>1
        for ii = 1:numel(varargin)/2
            ex.((varargin{2*ii-1})) = varargin{2*ii};
        end
    end

    % keyboard
    %% integrate the model forward
    for i = 1:numSteps
        if mod(i/numSteps, 0.05) == 0
            fprintf('%.2f%%\n', i/numSteps*100)
        end
        ex = ex.integrate;

        % save variables
        for j = 1:numel(variablesToSave)
            ExmirasRun.(variablesToSave{j})(i,:) = squeeze(ex.(variablesToSave{j})(1,1,:));
        end
        ExmirasRun.Zdr(i,:) = ex.Zhh(1,1,:)-ex.Zvv(1,1,:);

        if mod(i, 120) == 0
            fprintf('saving dsd at time step %d\n', i)
            ExmirasRun.N(end+1,:,:) = squeeze(ex.N(1,1,:,:));
            ExmirasRun.dNevap(end+1,:,:) = squeeze(ex.dNevap(1,1,:,:));
            ExmirasRun.dNcoal(end+1,:,:) = squeeze(ex.dNcoal(1,1,:,:));
            ExmirasRun.dNfallout(end+1,:,:) = squeeze(ex.dNfallout(1,1,:,:));
            % ExmirasRun.dN = squeeze(ex.dN(1,1,:,:
        end

        % stop integrating if the water vapor saturation is abov 98%
        if all(ex.qv(1,1,:)>=0.98)
            % remove all unused time steps
            for j = 1:numel(variablesToSave)
                ExmirasRun.(variablesToSave{j})(i+1:end,:) = [];
            end
            ExmirasRun.Zdr(i+1:end,:) = [];

            % end integration
            fprintf('Water vapor saturation reached at time step %d\n', i)
            break
        end 
    end

    % set up grids/save run
    ExmirasRun.Grids.tgrid = ((1:i))*ex.dt;
    ExmirasRun.Grids.zgrid = ex.zgrid;
    save([saveDir,ExmirasRun.ID,'.mat'], 'ExmirasRun')
    % load([saveDir,'S_ideal_300_0.30_40_0.50','.mat'], 'ExmirasRun')

    % % Plot the height-time profiles for the chosen variables
    % [TT, ZT] = meshgrid(ExmirasRun.Grids.tgrid, ExmirasRun.Grids.zgrid);

    % [~,dTdt] = gradient(ExmirasRun.T, 0.5); % convert to K/hr
    % ExmirasRun.dTdt = dTdt*3600;
    % ExmirasRun.RH = ExmirasRun.qv*100;
    % variablesToPlot = {'dTdt', 'RH', 'Zhh', 'Zdr', 'RR', 'rhohv', 'kdp'};
    % titles = {'CoolingRate', 'water vapor saturation', 'Z_{HH}', 'Z_{DR}', 'Rain Rate', 'Correlation Coefficient', 'Specific Differential Phase'};

    % plotTimeProfile(ExmirasRun.Grids.tgrid, ExmirasRun.Grids.zgrid, ['~/figures/publications/exmiras/',ExmirasRun.ID], ...
    %     cellfun(@(x) ExmirasRun.(x), variablesToPlot, 'uniform', 0), variablesToPlot)
    % cLabels = {'K', '', 'dBZ', 'dB', 'mm/hr', ' ', 'deg/km'};
    % clims = {[270,310], [0.1,1], [-2*10/6,50], [0,8], [0, 8], [0.99, 1], [0.5e-2, 0.5e-1]};
    % colormaps = {colortables('newblue'), flipud(winter(64)), colortables('RADAR32'), flipud(hsv(64)), flipud(spring(64)), flipud(turbo(64)), flipud(autumn(64))+[0,0,0.8]};
    % for i = 1:numel(variablesToPlot)
        % fig = figure("Units", "inches", "Position", [0,0,6.05,1.5]);

        % nuimagesc(TT, ZT, real(ExmirasRun.(variablesToPlot{i}))');
        % isc=imagesc(gca, ZT, TT, real(variablesToPlot{i}'));
        % isc.AlphaData = ones(size(isc.CData));
        % isc.AlphaData(isnan(real(variablesToPlot{iPlot}'))) = 0;
        % shading interp
        % xlabel('t [s]')
        % ylabel('z [m]')
        % title(titles{i})
        % % cb = colorbar;
        % setMeteoColormap(gca, variablesToPlot{i})
        % cb.Label.String = cLabels{i};
        % colormap(colormaps{i})
        % clim(clims{i})
        % shading flat
        % titles{i} = strrep(strrep(titles{i}, "{", ''), "}", '');
        % print2(fig, strrep(sprintf('~/figures/publications/exmiras/%s-%s.pdf', ExmirasRun.ID, variablesToPlot{i}), ' ', '-'))
    % end
    toc
end % end of runEXMIRAS

