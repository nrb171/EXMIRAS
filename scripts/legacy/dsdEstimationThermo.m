d = "/h/eol/nbarron/work/precip/spol-rhi/20220610";

files = dir(d+"/*Rhi1*.nc");
% f = "cfrad.20220610_044259.589_to_20220610_044516.729_SPOL_PrecipRhi2_RHI.nc"
% f = "cfrad.20220610_043616.597_to_20220610_043648.125_SPOL_PrecipRhi1_RHI.nc";;

for ii = 1:numel(files)
    f = files(ii).name;
    nci = ncinfo(d+"/"+f);

    range = ncread(d+"/"+f, 'range');
    azimuth = ncread(d+"/"+f, 'azimuth');
    elevation = ncread(d+"/"+f, 'elevation');


    % zhh and zdr flattened
    zhhf = ncread(d+"/"+f, 'DBZ');

    zdrf = ncread(d+"/"+f, 'ZDR_F');

    % reshape to 2D
    zhh = reshape(zhhf, [], numel(azimuth));
    zdr = reshape(zdrf, [], numel(azimuth));

    
    % zdr(zhh < 10) = NaN;
    % zhh(zhh < 10) = NaN;
    
    [~,rmesh] = meshgrid(elevation, range);
    zmesh = sind(elevation)'.*rmesh;
    
    figure
    subplot(2,1,1)
    pcolor(rmesh/1000, zmesh, zhh)
    shading flat
    ylim([0,12000])
    xlim([50, 150])
    colorbar
    
    setMeteoColormap(gca, 'Zhh')
    
    title(datestr(datetime(files(ii).name(7:21), 'InputFormat', 'yyyyMMdd_HHmmss')))
    
    % set(gca, 'YDir', 'normal')
    % print2
    
    
    % input('press enter to continue')
    % figure
    subplot(2,1,2)
    pcolor(rmesh/1000, zmesh, zdr)
    shading flat
    colorbar
    ylim([0,12000])
    xlim([50, 150])
    setMeteoColormap(gca, 'Zdr')
    % 20220610_040016
    
    % set(gca, 'YDir', 'normal')
    print2(gcf, sprintf('/h/eol/nbarron/figures/precip/spol-rhi/rhi-overview-%s.png', files(ii).name(7:21)))


    % input('press enter to continue')
end

ii=26
f = files(ii).name;
nci = ncinfo(d+"/"+f);

range = ncread(d+"/"+f, 'range');
azimuth = ncread(d+"/"+f, 'azimuth');
elevation = ncread(d+"/"+f, 'elevation');


% zhh and zdr flattened
zhhf = ncread(d+"/"+f, 'DBZ');

zdrf = ncread(d+"/"+f, 'ZDR_F');

% reshape to 2D
zhh = reshape(zhhf, [], numel(azimuth));
zdr = reshape(zdrf, [], numel(azimuth));



mask = ~isnan(zhh) & ~isnan(zdr);
zhhf = scatteredInterpolant(double(rmesh(mask)), double(zmesh(mask)), double(10.^(zhh(mask)/10)), 'natural', 'none');
zdrf = scatteredInterpolant(double(rmesh(mask)), double(zmesh(mask)), double(zdr(mask)), 'natural', 'none');

[zmeshp, rmeshp] = meshgrid(25:50:2025,65000:500:68000);
zhhp = 10*log10(zhhf(rmeshp, zmeshp));
zdrp = zdrf(rmeshp, zmeshp);

figure
hold on
colors = parula(size(zhhp,1));
for i = 1:size(zhhp,1)
    plot(zhhp(i,:), zmeshp(i,:), 'Color', colors(i,:))
end

print2


figure
hold on
% colors = parula(size(zhhp,1));
for i = 1:size(zhhp,1)
    plot(zdrp(i,:), zmeshp(i,:), 'Color', colors(i,:))
end

print2