classdef dsdAssimilation < handle
    properties
        bandName = 'S';
        
        zgrid = 25:50:2025;
        RHProfileObs = ones(1,41)*0.3;
        ZhhProfileObs = ones(1,41)*30;
        ZdrProfileObs = ones(1,41)*0.5;

        ZhhProfileSim 
        ZdrProfileSim 

        ZhhGrid = 20:5:75;
        DmGrid = 0.1:0.1:1.8;
        ZdrGrid = 0.3:0.3:3;

        RHGrid
        RHGridlores = 0.1:0.3:0.9;
        RHGridhires = 0.1:0.05:0.95;



        Ntop
        deta

        hires = true;
        minuteIdx= 10


        %% other classes used
        ra

        dNProfile

    end
    methods % dependent properties

        function deta = get.deta(obj)
            if obj.hires
                detat = load(sprintf('../data/LUTs/detaLUTs-hires-%1.0f.mat', obj.minuteIdx), 'deta');
                deta = detat.deta;
            else
                detat = load('../data/LUTs/detaLUTs-%1.0f.mat', 'deta');
                deta = detat.deta;
            end
        end

        function RHGrid = get.RHGrid(obj)
            if obj.hires
                RHGrid = obj.RHGridhires;
            else
                RHGrid = obj.RHGridlores;
            end
        end
    end

    methods % main methods
        
        function obj = getNtop(obj)

            % detaf = griddedInterpolant({zdrgrid, RHgrid, Dmgrid}, deta.(bandName), 'linear', 'extrap');

            %% Get Ntop for each combination of 
            obj.ex = exmiras;
            for zdr = obj.ZdrGrid
                for Dm = obj.DmGrid
                    
                    % Run the ideal simulation                  
                    %%! set up initial droplet size distribution
                    % fprintf('T0m = %1.2f, humidity = %1.2f, dBZStart = %d, zdr = %1.2f\n', T0m, humidity, dBZStart, zdr)
                    obj.ex.rngToggle = false;
                    %% set up initial droplet size distribution
                    obj.ex.xgrid = [50];
                    obj.ex.ygrid = [50];
                    obj.ex.zgrid = 25:50:2025;

                    sx = numel(obj.ex.xgrid);
                    sy = numel(obj.ex.ygrid);
                    sz = numel(obj.ex.zgrid);
                    Zhhi = obj.ZhhProfileObs(end) + zeros(sx, sy, sz)*0;
                    Zdri = obj.ZdrProfileObs(end) + rand(sx, sy, sz)*0;
                    Zhhi(:,:,1:size(Zhhi,3)-1) = -inf;
                    Zdri(:,:,1:size(Zhhi,3)-1) = 0;

                    obj.ex = obj.ex.initFromLambdaName(obj.bandName);
                    obj.ex = obj.ex.initFromDm(Zhhi(end), Zdri(end), Dm);
                    
                    Ntop3(...
                        obj.ZdrGrid ==zdr, ...
                        :, ...
                        obj.DmGrid == Dm, ...
                    :) = repmat(squeeze(obj.ex.N(1,1,end,:)),[1,3])';

                end
                
            end
            obj.Ntop = Ntop3;
        end

        

        function obj = estimateRadarProfiles(obj)
            %% OBSOLETE, USE estimateSingleRadarProfile INSTEAD AND THE DA METHODS BELOW

            %% initialize the interpolants
 
            Ntopf = griddedInterpolant({obj.ZdrGrid, obj.RHGrid, obj.DmGrid}, obj.Ntop, 'linear', 'linear');
            detaf = griddedInterpolant({obj.ZdrGrid, obj.RHGrid, obj.DmGrid}, obj.deta.(obj.bandName), 'linear', 'linear');



            np = NaN(numel(obj.zgrid), numel(obj.DmGrid), 250);


            Ntop3s = [];
            %% calculate the difference kernel at each height, RH, Dm (RH doesn't affect the Ntop)
            for kk = 1:numel(obj.zgrid)
                for ii = 1:numel(obj.DmGrid)
                    np(kk,ii,:) = detaf(obj.ZdrProfileObs(end), obj.RHProfileObs(kk), obj.DmGrid(ii));
                    Ntop3s(kk,ii,:) = Ntopf(obj.ZdrProfileObs(end), obj.RHProfileObs(1), obj.DmGrid(ii));
                end
            end

            

            % integrate the difference kernel to get the DSD at each height
            dnint = flipud(cumtrapz(obj.zgrid, np, 1));
            obj.dNProfile = Ntop3s + dnint.*max(Ntop3s, [], 3);

            for jj = 1:numel(obj.DmGrid)
                % calculate the radar observables from the DSD at each height
                for kk = 1:numel(obj.zgrid)
                    % print2
                    obj.ex.N(1,1,kk,:) = squeeze(obj.dNProfile(kk,jj,:))';
                    
                    obj.ZhhProfileSim(kk,jj) = obj.ex.Zhh(1,1,kk);
                    obj.ZdrProfileSim(kk,jj) = obj.ex.Zdr(1,1,kk);
                end
            end

            obj.ZhhProfileSim(imag(obj.ZhhProfileSim)>0) = -inf;
            obj.ZdrProfileSim(imag(obj.ZhhProfileSim)>0) = 0;
        end


        function varargout = estimateSingleRadarProfile(obj, Zhh, Zdr, Dm)
            % this is our backwards model that estimates the DSD profile.
            % in: Zhh in dBZ, Zdr in dB, Dm in mm for the top of the domain
            % out: dN, Zhhp, Zdrp, (dnint)

            %% load the Ntop LUT and interpolate to get the Ntop
            f=load('../data/LUTs/NtopLUTs.mat', 'NtopLUT', 'ZhhGrid', 'ZdrGrid', 'DmGrid');
            Ntopf = griddedInterpolant({10.^(f.ZhhGrid./10), f.ZdrGrid, f.DmGrid}, f.NtopLUT.(obj.bandName), 'linear', 'linear');
            Ntop = Ntopf(10.^(Zhh./10), Zdr, Dm);

            %% initialize the difference kernel interpolant
            detaf = griddedInterpolant({obj.ZdrGrid, obj.RHGrid, obj.DmGrid}, obj.deta.(obj.bandName), 'linear', 'linear');

            %% calculate the difference kernel (detakk) at each height, RH
            detakk = squeeze(detaf(...
                Zdr.*ones(size(obj.RHProfileObs)), ...
                obj.RHProfileObs, ...
                Dm.*ones(size(obj.RHProfileObs)) ...
            ));

            %% repeat the ntop for each height
            Ntop3s = squeeze(repmat(Ntop, [numel(obj.zgrid), 1]));

            % integrate the difference kernel to get the DSD at each height
            dnint = flipud(cumtrapz(obj.zgrid, detakk.*max(Ntop3s, [], 2), 1));
            dN = Ntop3s + dnint;
            dN(dN<0) = 0; % make sure we don't have negative DSDs (happens in rare cases)

            %% calculate the radar observables from the DSD at each height
            for kk = 1:numel(obj.zgrid)
                Zhhp(kk) = obj.ra.calcZhh(dN(kk,:));
                Zdrp(kk) = obj.ra.calcZdr(dN(kk,:));
            end

            %% return the requested outputs
            if nargout == 3
                varargout = {dN, Zhhp, Zdrp};
            elseif nargout == 4
                varargout = {dN, Zhhp, Zdrp, dnint};
            end

        end

         %% initialize N based on reflectivity and Zdr
        function [N,N0, mu, gamma] = getNFromZhhZdr(da, Z, Zdr)
            % in: Z in dBZ, Zdr in dBz
            % out: N in m^-3 mm^-1

            inds = find(~isinf(Z));
            mu = NaN(size(Z));
            gamma = NaN(size(Z));
            
            for i = 1:length(inds)
                % only perform search when there is data
                ind = inds(i);
                fun = @(x) calcErr(da.ra, da.ra.dpp, Z(ind), Zdr(ind), initN(da.ra, da.ra.dpp, Z(ind), x(1), x(2)), x(1), x(2));
                
                % initial guesses
                gammai = max(min(6./(Zdr),12), 0); % initial guess for gamma
                mui = -0.0201*gammai.^2 + 0.902*gammai - 1.78;% initial guess for mu
                
                % search for the optimal gamma and mu
                [x,~] =  fminsearchbnd(fun, [gammai, mui], [0, -2], [20, 15]);

                N(inds(i), :) = initN(da.ra, da.ra.dpp, Z(ind), x(1), x(2));
                mu(inds(i)) = x(2);
                gamma(inds(i)) = x(1);
                N0(inds(i)) = getN0(da.ra, da.ra.dpp, Z, gamma, mu);
            end
            function N0 = getN0(radar, dpp, dBZi, gamma, mu)
                % various preprocessing steps to calculate dual-pol variables
                mu = mu(:);
                gamma = gamma(:);
                N1 = (radar.D.^(mu) .* exp(-gamma.*radar.D));
                Zhh = radar.calcZhh(N1, dpp);
                N0 = 10.^(dBZi(:)/10)./10.^(Zhh(:)/10);
                N0 = reshape(N0, size(radar.mu));
            end

            function [N] = initN(radar, dpp, dBZi, gamma, mu)
                % various preprocessing steps to calculate dual-pol variables
                N1 = (radar.D.^(mu) .* exp(-gamma*radar.D));
                Zhh = 10^(calcZhh(radar, N1, dpp)/10);
                N0 = 10^(dBZi/10)./(Zhh);
                N = N0(:) .* (radar.D.^(mu) .* exp(-gamma*radar.D));
            end
            function err = calcErr(radar, dpp, dBZi, Zdri, N, gamma, mu)

                Zhh2 = radar.calcZhh(N, dpp);
                Zvv2 = radar.calcZvv(N, dpp);
                dZdr = Zhh2 - Zvv2 - Zdri;
                err = abs(dZdr);
                
            end
        end


        function [N,N0, mu, gamma] = getNFromZhhZdrDm(da,Z, Zdr, Dm)
            % in: Z in dBZ, Zdr in dBz
            % out: N in m^-3 mm^-1

            if any(isnan([Z, Zdr, Dm]))
                N = zeros(size(da.ra.D));
                N0 = NaN;
                mu = NaN;
                gamma = NaN;
                return
            end

            inds = find(~isinf(Z));
            [ix, iy, iz] = ind2sub(size(Z), inds);
            mu = NaN(size(Z));
            gamma = NaN(size(Z));
            
            fun = @(x) calcErr(da.ra, da.ra.dpp, Z, Zdr, Dm, initN(da.ra, da.ra.dpp, Z, x(1), x(2)), x(1), x(2));
            if ~da.ra.rngToggle
                gammai = max(min(6./(Zdr),12) + (2*rand()-1), 0); % initial guess for gamma
                % mui = -0.016*gammai.^2 + 1.213*gammai - 1.957 + min(max((8*rand()-4), -2), 15);% initial guess for mu
                mui = -0.0201*gammai.^2 + 0.902*gammai - 1.78 + min(max((2*rand()-1), -2), 15);% initial guess for mu
                [x,~] =  fminsearchbnd(fun, [gammai, mui], [0, -2], [20, 15]);
            else
                % ensemble approach to test convergence around above solution
                for ii = 1:20
                    gammai = randi(15);
                    mui = randi(12)-2;
                    [xt,~] =  fminsearchbnd(fun, [gammai, mui], [0, -2], [20, 15]);
                    x(ii, :) = xt;
                end
                x = mode(x, 1);
            end

            mu = x(2);
            gamma = x(1);
            N0 = getN0(da.ra, da.ra.dpp, Z, gamma, mu);
            N = N0.*da.ra.D.^mu.*exp(-gamma.*da.ra.D);

            function N0 = getN0(radar, dpp, dBZi, gamma, mu)
                % find the initial N0 based on Zhh, gamma, and mu
                mu = mu(:);
                gamma = gamma(:);
                N1 = (radar.D.^(mu) .* exp(-gamma.*radar.D));
                Zhh = calcZhh(radar, N1);
                N0 = 10.^(dBZi(:)/10)./10.^(Zhh(:)/10);
            end

            function [N] = initN(radar, dpp, dBZi, gamma, mu)
                % find the initial N based on Zhh, gamma, and mu
                N1 = (radar.D.^(mu) .* exp(-gamma*radar.D));
                Zhh = 10.^(radar.calcZhh(N1)/10);
                N0 = 10^(dBZi/10)./(Zhh);
                N = N0(:) .* (radar.D.^(mu) .* exp(-gamma*radar.D));
            end
            function err = calcErr(radar, dpp, dBZi, Zdri, Dm, N, gamma, mu)
                % calculate the error between the observed and simulated Zhh and Zdr
                Zhh2 = radar.calcZhh(N);
                Zvv2 = radar.calcZvv(N);
                dZdr = Zhh2 - Zvv2 - Zdri;
                [~,dmind]=max(N);
                Dm2 = radar.D(dmind);
                err = abs(abs(dZdr) + 10*abs((Dm2 - Dm)));
                
            end
        
        end
        function [N] = getNFromN0MuLambda(obj, N0, mu, lambda)
            % in: N0 in m^-3 mm^-1, mu, and lambda mm^-1 
            % out: N in m^-3 mm^-1
            
            N = N0(:) .* (obj.ra.D.^(mu(:))) .* exp(-lambda(:).*obj.ra.D);
            N(isnan(N)) = 0;
        end
    end

    methods 
        function bakeDSDs(obj)
            % function that calculates the DSD at the top of the rainshaft for many combinations of Zhh, Zdr, Dm
            % and saves them in a LUT. These can then be used later to interpolate the DSD at the top of the rainshaft for quick access.

            obj.ex = exmiras;

            NtopLUT = struct(...
                'S', NaN(numel(obj.ZhhGrid), numel(obj.ZdrGrid), numel(obj.DmGrid), 250), ...
                'C', NaN(numel(obj.ZhhGrid), numel(obj.ZdrGrid), numel(obj.DmGrid), 250), ...
                'X', NaN(numel(obj.ZhhGrid), numel(obj.ZdrGrid), numel(obj.DmGrid), 250) ...
            );
            progress = 0;
            total = numel(obj.ZhhGrid) * numel(obj.ZdrGrid) * numel(obj.DmGrid);
            for bandName = ["S", "C", "X"]             
                for Zhh = obj.ZhhGrid
                    for Dm = obj.DmGrid
                        for Zdr = obj.ZdrGrid
                            obj.ex = exmiras;

                            if mod(progress/total*100, 10) == 0
                                fprintf('Progress: %d/%d\n', progress, total)
                            end
            
                            % Run the ideal simulation                  
                            %%! set up initial droplet size distribution
                            % fprintf('T0m = %1.2f, humidity = %1.2f, dBZStart = %d, zdr = %1.2f\n', T0m, humidity, dBZStart, zdr)
                            obj.ex.rngToggle = false;
                            %% set up initial droplet size distribution
                            obj.ex.xgrid = [50];
                            obj.ex.ygrid = [50];
                            obj.ex.zgrid = 25:50:2025;

                            sx = numel(obj.ex.xgrid);
                            sy = numel(obj.ex.ygrid);
                            sz = numel(obj.ex.zgrid);
                            Zhhi = Zhh + zeros(sx, sy, sz)*0;
                            Zdri = Zdr + rand(sx, sy, sz)*0;
                            Zhhi(:,:,1:size(Zhhi,3)-1) = -inf;
                            Zdri(:,:,1:size(Zhhi,3)-1) = 0;

                            obj.ex = obj.ex.initFromLambdaName(bandName);
                            obj.ex = obj.ex.initFromDm(Zhhi(end), Zdri(end), Dm);
                            
                            NtopLUT.(bandName)(...
                                obj.ZhhGrid ==Zhh, ...
                                obj.ZdrGrid == Zdr, ...
                                obj.DmGrid == Dm, ...
                                : ...
                            ) = squeeze(obj.ex.N(1,1,end,:));

                            progress = progress + 1;
                        end
                    end
                end
            end

            ZhhGrid = obj.ZhhGrid;
            ZdrGrid = obj.ZdrGrid;
            DmGrid = obj.DmGrid;
            save('../data/LUTs/NtopLUTs.mat', 'NtopLUT', 'ZhhGrid', 'ZdrGrid', 'DmGrid')

        end
    end

    %% data assimilation methods
    methods
        function varargout=profileOptimizer(obj, ZhhProfileObs, ZdrProfileObs)
            % function that optimizes a profile of simulated Zhh and Zdr to find the best fitting DSD.
            % inputs: ZhhProfileObs, ZdrProfileObs: observed profiles of Zhh and Zdr to fit
            % outputs: ZhhOpt, ZdrOpt, DmOpt: optimized values of Zhh, Zdr, and Dm
            %          fv: value of the objective function at the optimum 
            

            
            options = optimset();
            options.TolX = 1e-2;
            options.TolFun = 1e-2;
        
            %% run the assimilation optimization routine
            % initial guess: median of the observed profiles, and 0.7 mm for D
            try
                [x,fv]=fminsearchbnd(@errorFunc, ...
                    [prctile(ZhhProfileObs(:),50), prctile(ZdrProfileObs(:),50), 0.7], ...
                    [prctile(ZhhProfileObs(:),5), prctile(ZdrProfileObs(:),15), 0.4], ...
                    [prctile(ZhhProfileObs(:),95), prctile(ZdrProfileObs(:),85), 1.6], ...
                    options);
            catch
                warning('Optimization failed, returning NaNs')
                switch nargout
                    case 3
                        varargout = {NaN, NaN, NaN};
                    case 4
                        varargout = {NaN, NaN, NaN, NaN};
                end
                return
            end 

            %% extract the optimized and return
            ZhhOpt = x(1); % dBZ
            ZdrOpt = x(2); % dB
            DmOpt = x(3); % mm
            switch nargout
                case 3
                    varargout = {ZhhOpt, ZdrOpt, DmOpt};
                case 4
                    varargout = {ZhhOpt, ZdrOpt, DmOpt, fv};
            end

            function errorReturn=errorFunc(x)
                % calculates the error between the observed and simulated profiles for a given set of Zhh, Zdr, and Dm
                % x(1) = Zhh, x(2) = Zdr, x(3) = Dm
                Zhh = x(1); % dBZ
                Zdr = x(2); % dB
                Dm = x(3); % mm
                % function to optimize the profile of the DSD
                [dN, Zhhp, Zdrp] = obj.estimateSingleRadarProfile(Zhh, Zdr, Dm);
                
                ZhhProfileObs = 10.^(ZhhProfileObs./10);
                Zhhp = 10.^(Zhhp./10);
                errorRaw = (Zhhp(:) - ZhhProfileObs).*(ZhhProfileObs./max(ZhhProfileObs(:))).^2;
                errorMode = mean(errorRaw, 2, 'omitnan');
                errorRawZdr = Zdrp(:) - ZdrProfileObs;
                errorModeZdr = mean(errorRawZdr, 2, 'omitnan');
                errorReturn = abs(mean(errorMode, "all", 'omitnan')) + ...
                    abs(mean(errorModeZdr, "all", 'omitnan'));
                
                if isnan(errorReturn)
                    error('not enough data to constrain the optimization');
                end

            end
        end

        function varargout=pointOptimizer(obj, ZhhPointObs, ZdrPointObs)
           %% uses the Zhang et al. (2001) method to get an initial guess for N0, mu, and gamma
           % this is the same as getNFromZhhZdr, but here for clarity
            [N, N0, mu, gamma] = obj.getNFromZhhZdrDm(ZhhPointObs, ZdrPointObs, 0.7);
            if nargout == 4
                varargout = {N, N0, mu, gamma};
            elseif nargout == 1
                varargout = {N};
            end
        end
    end

    methods 
        function da = dsdAssimilation(lambda)
            % loads an instance of the radar class for the specified wavelength
            da.ra = radar(lambda);

        end
    end
end
