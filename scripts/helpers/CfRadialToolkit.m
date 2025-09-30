%
% CfRadialToolkit - A collection of static methods to work with CfRadial files
classdef CfRadialToolkit
    methods (Static)
        function [rmesh, zmesh,varout] = regridVariable(range, azimuth, elevation, varin)
            % quick tool to regrid from 1D CfRadial conventions to 2D meshgrid, 
            % useful for plotting

            % generate meshgrid for plotting
            [~,rmesh] = meshgrid(elevation, range);
            zmesh = sind(elevation)'.*rmesh;
            rmesh = cosd(elevation)'.*rmesh;
            
            varout = reshape(varin, [], numel(azimuth));
        end

        function [] = plotCfRadial(files, varargin)
            %% Inputs
            % files should be a struct array as returned by dir
            % 'VariableName', VariableString, VariableString is a string array

            p = inputParser;
            addParameter(p, 'VariableName', 'DBZ');
            addParameter(p, 'Units', 'inches');
            addParameter(p, 'Position', [0,0,3,1.5]);
            addParameter(p, 'MeteoColormapName', 'Zhh');
            addParameter(p, 'SaveDirectory', '/h/eol/nbarron/figures/precip/spol-rhi/');
            addParameter(p, 'ylim', []);
            addParameter(p, 'xlim', []);
            addParameter(p, 'Platform', 'SPOL');
            addParameter(p, 'Azimuth', NaN);
            parse(p, varargin{:});
            
            if numel(p.Results.VariableName) ~= numel(unique(p.Results.MeteoColormapName))
                error('Number of variables and number of colormaps must be the same')
            end

            for ii = 1:numel(files)
                for jj = 1:numel(p.Results.VariableName)
                    for azi = p.Results.Azimuth
                        variable = p.Results.VariableName(jj);
                        f = files(ii).name;
                        d = files(ii).folder;

                        % load geometry and variable
                        
                        % keyboard
                        range = ncread(fullfile(files(1).folder,files(1).name), 'range');
                        azimuth = ncread(fullfile(files(1).folder,files(1).name), 'azimuth');
                        elevation = ncread(fullfile(files(1).folder,files(1).name), 'elevation');
                        varf = ncread(fullfile(files(1).folder,files(1).name), variable);
                        if ~isnan(azi)
                            mask = abs(azimuth-azi)<5; % only look at azimuths near the specified azimuth
                        else
                            mask = true(size(azimuth)); % no azimuth filtering
                        end

                        [elevmesh,rangemesh] = meshgrid(elevation(mask), range);
                        xmesh = rangemesh.*cosd(elevmesh);
                        ymesh = rangemesh.*sind(elevmesh);
                        var = varf(:,mask);

                        % plot figure
                        fig=figure('Units', p.Results.Units, 'Position', p.Results.Position);

                       
                        pcolor(gca,xmesh/1000, ymesh/1000,  var, 'EdgeColor', 'none')
                 
                        try
                            setMeteoColormap(gca,p.Results.MeteoColormapName(jj))
                        catch ME
                            ME.message
                            warning('Could not set meteo colormap, using default colormap instead')
                        end

                        % set limits and print
                        if ~isempty(p.Results.ylim)
                            ylim(p.Results.ylim)
                        end
                        if ~isempty(p.Results.xlim)
                            xlim(p.Results.xlim)
                        end
                        % keyboard

                        if ~isnan(azi)
                            print2(fig, sprintf('%srhi-%s-%s-%1.0f.pdf', p.Results.SaveDirectory, files(ii).name(7:21), variable, azi), 'Quality', '-r200')
                        else
                            print2(fig, sprintf('%srhi-%s-%s-all.pdf', p.Results.SaveDirectory, files(ii).name(7:21), variable), 'Quality', '-r200')
                        end
                    end
                end
            end
        end
    end
end