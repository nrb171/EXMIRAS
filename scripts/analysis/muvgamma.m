ex = exmiras;

gamma = [];
mu = [];
zdr1 =[];
N0 = [];
colors = winter(100);
zdrs = linspace(0.2, 5, 100);
zhh=30;
for ii = 1:numel(zdrs)
    zdr = zdrs(ii);
    ex = ex.initFromReflectivity(zhh, zdr);
    gamma(end+1) = ex.gamma;
    mu(end+1) = ex.mu;
    zdr1(end+1) = zdr;
    N0(end+1) = ex.N00;
end



% figure
% gamma = @(zdr) 10./(zdr)
% plot(gamma(zdr1), zdr1, 'LineWidth', 2);
% xlim([0,10])
% print2

mask = zdrs < 4;
pf = polyfit(gamma(mask), mu(mask), 2);
mufun = @(gamma)-0.016*gamma.^2 + 1.213*gamma -1.957;

fig = figure("Units","inches","Position", [0,0,3.25,3]);
scatter(gamma, mu, [],colors, 'filled');
hold on;
plot(gamma, mufun(gamma), 'r', 'LineWidth', 2);
plot(gamma(mask), polyval(pf, gamma(mask)), 'k--', 'LineWidth', 2);
xlabel('\Lambda');
ylabel('\mu');
colorbar;
colormap(colors);
clim([0.2,4.5])
print2(fig, 'muvgamma.png');
N2=[];
for i = 1:numel(gamma)
    D = 0.2:0.01:8;
    N = ((D.^(mu(i)) .* exp(-gamma(i)*D)));
    N(isinf(N)) = 0; % remove infs
    N = N./max(N); % normalize
    if sum(N) == 0
        keyboard
    end
    N2(i,:) = N;
end
% plot(D, N, 'Color', colors(i,:), 'LineWidth', 2);

fig = figure("Units","inches","Position", [0,0,3.25,3]);
hold on

set(gca,"Position", [0.1,0.05, 0.65, 0.95]);
pcolor(D, zdrs, N2)
shading flat
ylabel('Z_{DR} (dB)');
xlabel('D (mm)');
xticks(gca, 0:2:7);
ax2=axes("Position", [0.75,0.05, 0.2, 0.95]);

plot(N0/max(N0),zdrs)
xlim([0.000000001, 1])
ylim([0.2, 5])
xticks([0.000000001,1])
xscale("log")
set(ax2,'YTickLabel', repmat("", 1, numel(ax2.YTick)));
print2(fig, 'N.png');

% fig = figure("Units","inches","Position", [0,0,3.25,3]);
% scatter(gamma, zdr1)
% ylabel('ZDR');
% xlabel('\Gamma');
% print2(fig, 'zdrmu.png');

fig = figure("Units","inches","Position", [0,0,3.25,3]);
scatter(gamma, N0)
ylabel('N_0');
xlabel('\Lambda');
yscale('log');
print2(fig, 'n0Gamma.png');

i=1;



% reflh = readmatrix('../data/LUTs/reflh.txt');
% reflv = readmatrix('../data/LUTs/reflv.txt');

% plot(D,  0.9951 + 0.02510*D - 0.03644*D.^2 + 0.005303*D.^3 - 0.0002492*D.^4);
