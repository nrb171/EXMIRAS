



%%% 1. process thermodynamics profile %%%
ua=readmatrix("/h/eol/nbarron/workshop/EXMIRAS/test/data/hsinchu-2022061006.csv");
p = ua(:,3); % hPa
T = ua(:,4); % C
dwp = ua(:,12); % C
RH = ua(:,5); % %
Z = ua(:,10); % m

figure('Units', 'inches', 'Position', [0,0,3,3])
skewt(p, T, dwp)
print2('~/figures/publications/exmiras/hsinchu-2022061006-skewt.pdf')

zgrid = 25:50:2025;
RHProfile = interp1(Z, RH, zgrid, "linear", "extrap")/100;

de = dsdEstimator;
de.ZhhProfileObs = ones(1,41)*30;
de.ZdrProfileObs = ones(1,41)*0.3;
de.RHProfileObs = ones(1,41)*0.5;
de.RHProfileObs = RHProfile;
de.bandName = 'C';
de = de.getNtop();
de=de.estimateRadarProfiles();


ex = exmiras;
ex = ex.initFromLambdaName('S');
ex.N = zeros(1,1,1,numel(ex.D));
Zhh = [];
Zdr = [];
ZhhCalculated = [];
colors = winter(numel(de.DmGrid));
colors2 = autumn(numel(de.DmGrid));
fig1 = figure('Units', 'inches', 'Position', [0,0,3,3]);
ax1 = axes(fig1);
hold on
fig2 = figure('Units', 'inches', 'Position', [0,0,3,3]);
ax2 = axes(fig2);
hold on

fig3 = figure('Units', 'inches', 'Position', [0,0,3,3]);
ax3 = axes(fig3);
hold on
for jj = 1:numel(de.DmGrid)
    hold on
    

    plot(ax1,de.ZhhProfileSim(:,jj), zgrid, 'Color', colors(jj,:))
    % hold on
    % plot(ZhhCalculated(:,jj), zgrid, 'Color', colors2(jj,:))
    % subplot(1,2,2)
    plot(ax2,de.ZdrProfileSim(:,jj), zgrid, 'Color', colors(jj,:))

    plot(ax3, de.ex.D, squeeze(sum(de.dNProfile(:,jj,:),1)), 'Color', colors(jj,:))
end

xlim(ax2, [round(min(de.ZdrProfileSim(:)),2)-0.1, round(max(de.ZdrProfileSim(:)),2)+0.1])
xlim(ax1, [floor(min(de.ZhhProfileSim(:)))-0.5, ceil(max(de.ZhhProfileSim(:)))+0.5])
print2(fig1, '~/figures/publications/exmiras/hsinchu-2022061006-Zhh-vs-height.pdf')

print2(fig2, '~/figures/publications/exmiras/hsinchu-2022061006-Zdr-vs-height.pdf')


xlim(ax3, [0,4])
ylim(ax3, [-1e5, 1e5])
print2(fig3, '~/figures/publications/exmiras/hsinchu-2022061006-N-vs-D.pdf')



figure
hold on

