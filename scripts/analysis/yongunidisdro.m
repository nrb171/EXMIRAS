d = "/h/eol/nbarron/work/precip/disdro/";

files = dir(d+"*.nc");

figure
hold on
for ii = 1:numel(files)
    f = files(ii).name;
    d = files(ii).folder;
    
    RR = ncread(fullfile(d,f), 'precip_rate');
    tt = ncread(fullfile(d,f), 'time'); 
    yyyy = str2double(f(14:17));
    mm = str2double(f(18:19));
    dd = str2double(f(20:21));
    dt = minutes(tt) + datetime(yyyy,mm,dd);
    % title(f)
    plot(dt, RR, '-k');
    % pause
end

xlim([datetime(2022, 6,16,12,0,0), datetime(2022, 6,16, 16,0,0)])
print2

d = "/h/eol/nbarron/work/precip/disdro/";
f = "DSD_Yonaguni_20220616_qc.nc";

RR = ncread(fullfile(d,f), 'precip_rate');
[~,im]=max(RR);
D = ncread(fullfile(d,f), 'particle_size');
ndp = ncread(fullfile(d,f), 'qc_number_detected_particles');
rs = ncread(fullfile(d,f), 'raw_spectrum');

tt = ncread(fullfile(d,f), 'time');

figure
plot(D, sum(rs(:,:,im),2))
xlim([0 8])
print2
