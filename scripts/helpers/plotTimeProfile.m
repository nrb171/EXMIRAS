function plotTimeProfile(TT, ZT, prefix, variablesToPlot, titles, cLabels, clims, colormaps)
    % TT = (1:numSteps)*ex.dt;
    % ZT = ex.zgrid;
    % prefix = prefix to save the plots.
    % variablesToPlot = cell array of ZT x TT variables, e.g., {T, qv, dBZhh, Zdr, RR, rhohv, kdp};
    % titles = cell array of titles for each plot, e.g., {'T', 'qv', 'dBZhh', 'Zdr', 'RR', 'rhohv', 'kdp'};
    % cLabels = cell array of colorbar labels for each plot, e.g., {'T', 'qv', 'dBZhh', 'Zdr', 'RR', 'rhohv', 'kdp'};
    % clims = cell array of color limits for each plot, e.g., {[-5,5], [0,1], [-2*10/6,50], [-2*10/6,50], [0, 1], [0, 1], [0, 1]};
    % colormaps = cell array of colormaps for each plot, e.g., {hot(2), hot(2), hot(2), hot(2), hot(2), hot(2), hot(2)};


    %    % vertical profile vs time plot

    for iPlot = 1:numel(variablesToPlot)
        fig = figure("Units", "inches", "Position", [0,0,5.75,1.5]);

        pcolor(TT, ZT, real(variablesToPlot{iPlot}'))
        ylabel('z [m]')
        title(titles{iPlot})
        cb = colorbar;
        cb.Label.String = cLabels{iPlot};
        colormap(colormaps{iPlot})
        clim(clims{iPlot})
        xlim([datetime(2015,6,13,2,0,0), datetime(2015,6,13,6,0,0)])
        shading flat
        titles{iPlot} = strrep(strrep(titles{iPlot}, "{", ''), "}", '');
        print2(fig, strrep(sprintf('%s_model_%s.png', prefix,titles{iPlot}), ' ', '_'), 'Quality', '-r300')
    end
end