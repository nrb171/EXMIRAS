saveDir = '/h/eol/nbarron/work/exmirasData';


% try a bunch of initial DSDs and see how the difference kernel changes
for zdr=1:0.5:3
    for i=1:50
        runEXMIRAS(300, 0.3, 40, zdr, "C", -6.5, [saveDir, 'diff-kernel-', num2str(i)], 'waterVaporInteraction', false)
    end
end

return


files = dir([saveDir, 'diff-kernel-*40_0.5*.mat']);
% files = dir([saveDir, 'diff-kernel-*.mat']);

De = logspace(log10(0.1),log10(8), 250+1);
D =  (De(1:end-1) + De(2:end))/2;
figure
hold on
zdr = [];
a = [];
b = [];
zhh = [];
xs = [];
dm = [];
dnhat = [];
for i=1:numel(files)
    
    f = load([files(i).folder, '/', files(i).name]);

    
    
    Ntop = squeeze(f.ExmirasRun.N(end,end,:));
    Nbot = squeeze(f.ExmirasRun.N(end,1,:));
    Ndiff = Ntop - Nbot;
    dNevap = sum(mean(squeeze(f.ExmirasRun.dNevap(end,:,:)),1),1);
    dNcoal = sum(mean(squeeze(f.ExmirasRun.dNcoal(end,:,:)),1),1);
    % dNfallout = sum(squeeze(f.ExmirasRun.dNfallout

    gamma = @(N,mu,lambda) N.*D.^mu.*exp(-lambda.*D);
    funcTop = @(N,mu,lambda) rms(gamma(N,mu,lambda)-(Ntop)');
    foptTop = fminsearchbnd(@(x) funcTop(x(1), x(2), x(3)), [2.7408e+04,4, 3.5], [1e3, -2, 0], [1e9,15, 20]);
    funcBot = @(N,mu,lambda) rms(gamma(N,mu,lambda)-(Nbot)');
    foptBot = fminsearchbnd(@(x) funcBot(x(1), x(2), x(3)), [2.7408e+04,4, 3.5], [1e3, -2, 0], [4e9,15, 20]);
    funcDiff = @(N,mu,lambda) rms(gamma(N,mu,lambda)-((Nbot-Ntop)'));
    fopt = fminsearchbnd(@(x) funcDiff(-x(1), x(2), x(3)), foptTop, [foptTop(1)-0.1*foptTop(1), foptTop(2)-0.1*foptTop(2), 0], [foptTop(1)+0.1*foptTop(1), foptTop(2)+0.1*foptTop(2), 20]);
    fopt = fminsearchbnd(@(x) funcDiff(-x(1), x(2), x(3)), foptTop, [1e2, -2, 0], [4e9,15, 20])

    % asd = @(x) rms(exp(x(1).*D) + x(2)*D + x(3) - ((Nbot-Ntop)./Ntop)')

    
    % dn = @(x) rms(x(1)./(x(2).*D) + x(3)*D + x(4) - ((Nbot-Ntop)./Ntop)')
    Ndiff = ((Nbot-Ntop)./Ntop)';
    Ndiff(end-2:end) = Ndiff(end-3);
    % D(D>2) = [];
    param = @(x) +x(1).*exp(x(2).*D).*D.^x(3) +x(4).*exp(x(5).*D).*D.^x(6)+ x(7)
    % param = @(x) x(1)./(x(2).*D).*D.^x(3) + x(4)*D + x(5)
    % param = @(x) -x(1).*exp(x(2).*D).*D.^x(3) 
    % dn = @(x) rms(param(x) - Ndiff);
    
    
    % f2=fminsearchbnd(dn, [0.3874, -0.3492, 1.1881, -0.6005, -0.2303, 0.2147, -0.3663], [0, -3, -3, -1, -3, -3, -1],[1, 0, 3, 0, 0, 3, 1] )

    % figure
    % hold on
    % plot(D, Ndiff)
    % plot(D, param(f2))

    % ylim([-1,1])
    % print2

    % xs(1:7,end+1) = f2;


    % % ylim([-1,1])
    % hold on
    % % x=
    % figure
    % hold on
    % plot(D, gamma(-fopt(1),fopt(2), fopt(3))./gamma(foptTop(1),foptTop(2), foptTop(3)))
    plot(D, (Nbot-Ntop)./Ntop)
    % ylim([-1,1])
    % print2
    % print2

    [~,atemp] = max(Ntop);
    a(end+1) = D(atemp);
    [~, btemp] = min((Nbot-Ntop)./Ntop);
    b(end+1) = D(btemp);

    dm(end+1) = D(atemp);
    dnhat(end+1,:) = (Nbot-Ntop)./Ntop;
    % zdr(end+1) = f.ExmirasRun.Zdr(end,end);
    % zhh(end+1) = f.ExmirasRun.Zhh(end,end);
   
    % plot(D, -Ndiff./Ntop)
  
end
ylim([-1,1])
print2

figure
scatter(a,b)
xlim([0,2])
ylim([0,2])
print2

figure
hold on
mu = [];
lambda = [];
N0 = [];
md = [];
for ii = 1:300
    ex = exmiras;
    ex = ex.initFromLambdaName('C');
    ex.rngToggle = false;
    ex = ex.initFromReflectivity(40, 0.5);
    mu(end+1) = ex.mu;
    lambda(end+1) = ex.gamma;
    N0(end+1) = ex.N00;    
    md(end+1) =D(find(squeeze(ex.N(1,1,end,:)) == max(squeeze(ex.N(1,1,end,:)))));
    % plot(D, squeeze(ex.N(1,1,:)))
    % print2
end


figure('Units', 'inches', 'Position', [0,0,3,3])
hold on
scatter(md, mu, 'DisplayName', '\mu')
pfmu = polyfit(md, mu, 3)
plot(sort(md), polyval(pfmu, sort(md)), '-k', 'HandleVisibility', 'off')
scatter(md, lambda, 'DisplayName', '\lambda')
pflambda = polyfit(md, lambda, 3)
plot(sort(md), polyval(pflambda, sort(md)), '-k','HandleVisibility', 'off')
yyaxis right
scatter(md, N0, 'DisplayName', 'N_0')
pfN0 = polyfit(md, log10(N0), 3)
plot(sort(md), 10.^(polyval(pfN0, sort(md))), '-k','HandleVisibility', 'off')
% ylabel('N_0 [m^{-3} mm^{-1}]
set(gca, 'YScale', 'log')
xlabel('D_m [mm]')
legend
print2('~/figures/publications/exmiras/mu-lambda-N0-vs-Dm.pdf')
% func = @(N,mu,lambda) rms(-N.*D.^mu.*exp(-lambda.*D)-(Nbot-Ntop)');
% figure

% foptTop(1), foptBot(1), fopt(1)
% foptTop(2), foptBot(2), fopt(2)
% foptTop(3)-fopt(3)

% plot(D, -fopt(1)*D.^fopt(2).*exp

[~,ia]=unique(dm);
dnhatmap=interp1(dm(ia),dnhat(ia,:), 0.2:0.1:1.1, 'linear', 'extrap');
figure
pcolor(D, 0.2:0.1:1.1,dnhatmap)
clim([-1,1])
colorbar
shading flat
print2

figure
hold on
for ii = 1:size(dnhatmap,1)
    plot(D, dnhatmap(ii,:), 'DisplayName', sprintf('D_m = %1.1f mm', 0.1*ii+0.1))
end
ylim([-1,1])
legend
print2

dBZStart = 40;
zdr = 0.5;
ex = exmiras;
% ex.rngToggle = false;
%% set up initial droplet size distribution
ex.xgrid = [50];
ex.ygrid = [50];
ex.zgrid = 50:50:2050;
% ex = ex.initFromLambdaName('C');
sx = numel(ex.xgrid);
sy = numel(ex.ygrid);
sz = numel(ex.zgrid);
dBZi = dBZStart + rand(sx, sy, sz)*0;
Zdri = zdr + rand(sx, sy, sz)*0;
dBZi(:,:,1:size(dBZi,3)-1) = -inf;
Zdri(:,:,1:size(dBZi,3)-1) = 0;

ex = ex.initFromLambdaName('C');
ex = ex.initFromDm(40,0.5);

ex.NP = ex.N;

dms = 0.2:0.1:1.3;
Zhhmap = [];
zdrmap = [];
figure('Units', 'inches', 'Position', [0,0,3,3])
hold on
colors = winter(numel(dms));
Ni = [];
for i=1:numel(dms)
    ex = ex.initFromDm(40, 0.5,dms(i));
    Ntop = ex.N00;
    mu = ex.mu;
    lambda = ex.gamma;
    ex.Zdr(end)
    ex.N(1,1,41,:) = Ntop.*D.^mu.*exp(-lambda.*D);
    Ni(end+1) = interp1(D, squeeze(ex.N(1,1,41,:)), dms(i));
    meanD = trapz(D, Ntop.*D.^mu.*exp(-lambda.*D).*D)./trapz(D, Ntop.*D.^mu.*exp(-lambda.*D))
    plot(D, squeeze(ex.N(1,1,41,:)), "Color", colors(i,:), 'DisplayName', sprintf('D_m = %1.1f mm', dms(i)))
end
xlim([0,4])
set(gca, 'YScale', 'log')
plot(dms, Ni, '-k')
print2
print2('~/figures/publications/exmiras/N-vs-Dm.pdf')



figure('Units', 'inches', 'Position', [0,0,3,3])
% dnhatmap(:,end-3:end) = 
hold on
for ii = 1:size(Zhhmap,2)
    plot( D, dnhatmap(ii,:), 'DisplayName', sprintf('D_m = %1.1f mm', dms(ii)), 'Color', colors(ii,:))
end
% pcolor(dms, D, dnhatmap)
% shading flat
% set(gca, 'YScale', 'log')
% colorbar 
ylim([-1,1]) 
print2('~/figures/publications/exmiras/dnhat-vs-Dm.pdf')


figure('Units', 'inches', 'Position', [0,0,3,3])
hold on
for ii = 1:size(Zhhmap,2)
    plot(Zhhmap(:,ii), ex.zgrid,  'DisplayName', sprintf('D_m = %1.1f mm', dms(ii)), 'Color', colors(ii,:))
end
% pcolor(dms, ex.zgrid, Zhhmap)
% shading flat
% setMeteoColormap(gca, 'Zhh')
% legend
print2('~/figures/publications/exmiras/Zhh-vs-Dm.pdf')

figure('Units', 'inches', 'Position', [0,0,3,3])
hold on
for ii = 1:size(dnhatmap,2)
    plot(zdrmap(:,ii), ex.zgrid,  'DisplayName', sprintf('D_m = %1.1f mm', dms(ii)), 'Color', colors(ii,:))
end
% pcolor(dms, ex.zgrid, zdrmap)
% shading flat
% setMeteoColormap(gca, 'Zdr')
% legend
print2('~/figures/publications/exmiras/Zdr-vs-Dm.pdf')

