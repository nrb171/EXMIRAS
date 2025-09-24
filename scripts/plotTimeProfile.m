function plotTimeProfile(TT, ZT, prefix, variablesToPlot, titles)
    % TT = (1:numSteps)*ex.dt;
    % ZT = ex.zgrid;
    % prefix = prefix to save the plots.
    % variablesToPlot = cell array of ZT x TT variables, e.g., {T, qv, dBZhh, Zdr, RR, rhohv, kdp};
    % titles = cell array of titles for each plot, e.g., {'T', 'qv', 'dBZhh', 'Zdr', 'RR', 'rhohv', 'kdp'};
    % cLabels = cell array of colorbar labels for each plot, e.g., {'T', 'qv', 'dBZhh', 'Zdr', 'RR', 'rhohv', 'kdp'};
    % clims = cell array of color limits for each plot, e.g., {[-5,5], [0,1], [-2*10/6,50], [-2*10/6,50], [0, 1], [0, 1], [0, 1]};
    % colormaps = cell array of colormaps for each plot, e.g., {hot(2), hot(2), hot(2), hot(2), hot(2), hot(2), hot(2)};
% keyboard

    %    % vertical profile vs time plot
    ZT = interp1(0:numel(ZT)-1, ZT, -0.5:1:numel(ZT)-0.5, 'linear', 'extrap');
    TT = interp1(0:numel(TT)-1, TT, -0.5:1:numel(TT)-0.5, 'linear', 'extrap');
    for iPlot = 1:numel(variablesToPlot)
        fig = figure("Units", "inches", "Position", [0,0,6.05,1.5]);
        % nuimagesc(gca,TT, ZT, real(variablesToPlot{iPlot}'));
        % keyboard
        
        isc=imagesc(gca, ZT, TT, real(variablesToPlot{iPlot}'));
        isc.AlphaData = ones(size(isc.CData));
        isc.AlphaData(isnan(real(variablesToPlot{iPlot}'))) = 0;
        set(gca, 'YDir', 'normal');
        % ylabel('z [m]')
        % title(titles{iPlot})
        % cb = colorbar;
        % cb.Label.String = cLabels{iPlot};
        % colormap(colormaps{iPlot})
        % clim(clims{iPlot})
        setMeteoColormap(gca, titles{iPlot})
        % xlim([datetime(2015,6,13,2,0,0), datetime(2015,6,13,6,0,0)])
        % xlim([-2995, 11406])
        xlim([0,4500])
        xticks(0:1000:4000)
        ylim([0, 2500])
        shading flat
        titles{iPlot} = strrep(strrep(titles{iPlot}, "{", ''), "}", '');
        print2(fig, strrep(sprintf('%s_model_%s.pdf', prefix,titles{iPlot}), ' ', '_'), 'Quality', '-r300')
    end
end