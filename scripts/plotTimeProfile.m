function plotTimeProfile(TT, ZT, idName, variablesToPlot, titles, cLabels, clims, colormaps)
    % vertical profile vs time plot

    for iPlot = 1:numel(variablesToPlot)
        fig = figure("Units", "inches", "Position", [0,0,5.75,1.5]);

        pcolor(TT, ZT, real(variablesToPlot{iPlot}'))
        ylabel('z [m]')
        % t = strsplit(titles{i}, '_');
        title(titles{iPlot})
        cb = colorbar;
        cb.Label.String = cLabels{iPlot};
        colormap(colormaps{iPlot})
        clim(clims{iPlot})
        xlim([datetime(2015,6,13,2,0,0), datetime(2015,6,13,6,0,0)])
        shading flat
        titles{iPlot} = strrep(strrep(titles{iPlot}, "{", ''), "}", '');
        print2(fig, strrep(sprintf('~/figures/evap/%s_model_%s.png', idName,titles{iPlot}), ' ', '_'), 'Quality', '-r300')
    end
end