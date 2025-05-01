for T0m = 285:5:310
    for humidity = 0.1:0.1:0.9
        for dBZStart = 30:5:50
                    
            %%! set up initial droplet size distribution
            fprintf('T0m = %d, humidity = %.1f\n', T0m, humidity)
            ex = exmiras;
            %% set up initial droplet size distribution
            ex.xgrid = [50];
            ex.ygrid = [50];
            ex.zgrid = 50:50:2050;
            bandName = "S"

            sx = numel(ex.xgrid);
            sy = numel(ex.ygrid);
            sz = numel(ex.zgrid);
            zdr = 0.5
            dBZi = dBZStart + rand(sx, sy, sz)*0;
            Zdri = zdr + rand(sx, sy, sz)*0;
            dBZi(:,:,1:size(dBZi,3)-1) = -inf;
            Zdri(:,:,1:size(dBZi,3)-1) = 0;

            ex = ex.initFromLambdaName("S");
            ex = ex.initFromReflectivity(dBZi, Zdri);

            ex.NP = ex.N;

            %%! set up initial state variables
            temp = @(z) -6.5*z/1000 + T0m;
            pres = @(z) 1000*exp(-z/1000/8.4);
            es = @(T) 2.53e9*exp(-5420./(T)); %hPa
            [ym, xm, zm] = meshgrid(ex.ygrid, ex.xgrid, ex.zgrid);
            ex.T = temp(zm);
            ex.p = pres(zm);
            ex.pv = es(ex.T)*humidity;

            ex.u = zeros(size(ex.N, [1:3]))
            ex.v = zeros(size(ex.N, [1:3]))
            ex.w = zeros(size(ex.N, [1:3]))


            plot(ex.D, squeeze(ex.N(1,1,end,:)) )
            print2(gcf, './.temp.png')


            ex.dt = 0.5;
            ex.st = 3*3600./ex.dt; 

            fig = figure(units="inches", position=[0,0,3.75,3]);
            plot(ex.D, (squeeze(ex.dNcoal(1,1,end,:))) )
            hold on
            plot(ex.D, (squeeze(ex.dNevap(1,1,end,:)) ))
            yyaxis right
            plot(ex.D, squeeze(ex.dNfallout(1,1,end,:)))

            % plot(ex.D, (squeeze(ex.dNcoal(1,1,end,:)) )+squeeze(ex.dNevap(1,1,end,:)) )
            yyaxis left
            xscale("log")
            xlabel('D [mm]')
            ylabel('dN [m^{-3} mm^{-1}]')
            legend('coal', 'evap', 'fallout')
            print2(fig, './.temp.png')



            %% set storm time
            numSteps = ex.st+0;
            variablesToSave = {'T', 'p', 'pv', 'qv', 'Zhh', 'Zvv', 'RR', 'rhohv', 'kdp', 'Zdr', 'theta'};
            ExmirasRun = struct();
            ExmirasRun.initVariables.p = pres(zm);
            ExmirasRun.initVariables.T = temp(zm);
            ExmirasRun.initVariables.pv = es(ex.T)*humidity;
            ExmirasRun.ID = sprintf('%s_ideal_%d_%1.2f_%2.0d', bandName, T0m, humidity,dBZStart);
            for i = 1:numel(variablesToSave)
                ExmirasRun.(variablesToSave{i}) = NaN(numSteps, sz);
            end

            %% integrate the model forward
            for i = 1:numSteps
                if mod(i/numSteps, 0.05) == 0
                    fprintf('%.2f%%\n', i/numSteps*100)
                end
                ex = ex.integrate;

                % save variables
                for j = 1:numel(variablesToSave)
                    ExmirasRun.(variablesToSave{j})(i,:) = squeeze(ex.(variablesToSave{j})(1,1,:));
                end
                ExmirasRun.Zdr(i,:) = squeeze(10*log10(ex.Zhh./ex.Zvv));

                % stop integrating if the water vapor saturation is reached above 99%
                if all(ex.qv(1,1,:)>=0.95)
                    break
                end 
            end

            % set up grids/save run
            ExmirasRun.Grids.tgrid = ((1:i))*ex.dt;
            ExmirasRun.Grids.zgrid = ex.zgrid;
            save(['./',ExmirasRun.ID,'.mat'], 'ExmirasRun')

            % plotting
            [TT, ZT] = meshgrid(ExmirasRun.Grids.tgrid, ExmirasRun.Grids.zgrid);

            try
                variablesToPlot = {'T', 'qv', 'Zhh', 'Zdr', 'RR', 'rhohv', 'kdp'};
                titles = {'temperature perturbation', 'water vapor saturation', 'Z_{HH}', 'Z_{DR}', 'Rain Rate', 'Correlation Coefficient', 'Specific Differential Phase'};
                cLabels = {'K', '', 'dBZ', 'dB', 'mm/hr', ' ', 'deg/km'};
                clims = {[270,310], [0.1,1], [-2*10/6,50], [-8,8], [0, 8], [0.99, 1], [0.5e-2, 0.5e-1]};
                colormaps = {colortables('newblue'), flipud(winter(64)), colortables('RADAR32'), flipud(hsv(64)), flipud(spring(64)), flipud(turbo(64)), flipud(autumn(64))+[0,0,0.8]};
                for i = 1:numel(variablesToPlot)
                    fig = figure("Units", "inches", "Position", [0,0,5.75,1.5]);

                    pcolor(TT, ZT, real(ExmirasRun.(variablesToPlot{i}))');
                    shading interp
                    xlabel('t [s]')
                    ylabel('z [m]')
                    title(titles{i})
                    cb = colorbar;
                    cb.Label.String = cLabels{i};
                    colormap(colormaps{i})
                    clim(clims{i})
                    shading flat
                    titles{i} = strrep(strrep(titles{i}, "{", ''), "}", '');
                    print2(fig, strrep(sprintf('./%s-%s.png', ExmirasRun.ID, titles{i}), ' ', '-'))
                end
            catch ME
                fprintf('Error in plotting: %s\n', ME.message)
            end
        end
    end
end