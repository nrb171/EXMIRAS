% scan through all of the .cls files and find the mean RH between 0 and 2 km agl.

folder = "/h/eol/nbarron/work/pecan/upper-air/nssl";

files = dir(fullfile(folder, '**', 'NSSL2*.cls'));

for ii = 1:length(files)
    file = files(ii);
    tab = readmatrix(fullfile(file.folder, file.name), 'FileType', 'text', 'NumHeaderLines', 15);
    rh = tab(:,5);
    alt = tab(:,15) - tab(1,15); % altitude agl in m

    rh2km(ii) = mean(rh(alt <= 2000 & rh ~= 999));

end

rh2km(contains({files.name}, '0702'))

lb = 90
tab2 = table(string({files(rh2km > lb & rh2km<lb+10).name})',rh2km(rh2km > lb & rh2km<lb+10)', 'VariableNames', {'Filename', 'RH_0_2km'});
tab2
