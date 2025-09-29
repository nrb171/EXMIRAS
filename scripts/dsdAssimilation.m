classdef dsdAssimilation < handle
    properties
        bandName = 'S';
        D
        De
        Dw
        nBins = 250;
        dMax = 8; % mm

        
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
        detaInterpolant
        detaLUTs

        hires = true;
        minuteIdx= 7


        %% other classes used
        ra

        dNProfile

    end
    methods % dependent properties

        function deta = get.deta(obj)

            if isempty(obj.detaLUTs)

                if obj.hires
                    detat = load(sprintf('../data/LUTs/detaLUTs-hires-%1.0f.mat', obj.minuteIdx), 'deta');
                    deta = detat.deta;
                else
                    detat = load(sprintf('../data/LUTs/detaLUTs-%1.0f.mat', obj.minuteIdx), 'deta');
                    deta = detat.deta;
                end
                obj.detaLUTs = deta;
            else
                deta = obj.detaLUTs;
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
        
        % function obj = getNtop(obj)

        %     % detaf = griddedInterpolant({zdrgrid, RHgrid, Dmgrid}, deta.(bandName), 'linear', 'extrap');

        %     %% Get Ntop for each combination of 
        %     obj.ex = exmiras;
        %     for zdr = obj.ZdrGrid
        %         for Dm = obj.DmGrid
                    
        %             % Run the ideal simulation                  
        %             %%! set up initial droplet size distribution
        %             % fprintf('T0m = %1.2f, humidity = %1.2f, dBZStart = %d, zdr = %1.2f\n', T0m, humidity, dBZStart, zdr)
        %             obj.ex.rngToggle = false;
        %             %% set up initial droplet size distribution
        %             obj.ex.xgrid = [50];
        %             obj.ex.ygrid = [50];
        %             obj.ex.zgrid = 25:50:2025;

        %             sx = numel(obj.ex.xgrid);
        %             sy = numel(obj.ex.ygrid);
        %             sz = numel(obj.ex.zgrid);
        %             Zhhi = obj.ZhhProfileObs(end) + zeros(sx, sy, sz)*0;
        %             Zdri = obj.ZdrProfileObs(end) + rand(sx, sy, sz)*0;
        %             Zhhi(:,:,1:size(Zhhi,3)-1) = -inf;
        %             Zdri(:,:,1:size(Zhhi,3)-1) = 0;

        %             obj.ex = obj.ex.initFromLambdaName(obj.bandName);
        %             obj.ex = obj.ex.initFromDm(Zhhi(end), Zdri(end), Dm);
                    
        %             Ntop3(...
        %                 obj.ZdrGrid ==zdr, ...
        %                 :, ...
        %                 obj.DmGrid == Dm, ...
        %             :) = repmat(squeeze(obj.ex.N(1,1,end,:)),[1,3])';

        %         end
                
        %     end
        %     obj.Ntop = Ntop3;
        % end
        % function obj = estimateRadarProfiles(obj)
        %     %% OBSOLETE, USE estimateSingleRadarProfile INSTEAD AND THE DA METHODS BELOW

        %     %% initialize the interpolants
 
        %     Ntopf = griddedInterpolant({obj.ZdrGrid, obj.RHGrid, obj.DmGrid}, obj.Ntop, 'linear', 'linear');
        %     detaf = griddedInterpolant({obj.ZdrGrid, obj.RHGrid, obj.DmGrid}, obj.deta.(obj.bandName), 'linear', 'linear');



        %     np = NaN(numel(obj.zgrid), numel(obj.DmGrid), 250);


        %     Ntop3s = [];
        %     %% calculate the difference kernel at each height, RH, Dm (RH doesn't affect the Ntop)
        %     for kk = 1:numel(obj.zgrid)
        %         for ii = 1:numel(obj.DmGrid)
        %             np(kk,ii,:) = detaf(obj.ZdrProfileObs(end), obj.RHProfileObs(kk), obj.DmGrid(ii));
        %             Ntop3s(kk,ii,:) = Ntopf(obj.ZdrProfileObs(end), obj.RHProfileObs(1), obj.DmGrid(ii));
        %         end
        %     end

            

        %     % integrate the difference kernel to get the DSD at each height
        %     dnint = flipud(cumtrapz(obj.zgrid, np, 1));
        %     obj.dNProfile = Ntop3s + dnint.*max(Ntop3s, [], 3);

        %     for jj = 1:numel(obj.DmGrid)
        %         % calculate the radar observables from the DSD at each height
        %         for kk = 1:numel(obj.zgrid)
        %             % print2
        %             obj.ex.N(1,1,kk,:) = squeeze(obj.dNProfile(kk,jj,:))';
                    
        %             obj.ZhhProfileSim(kk,jj) = obj.ex.Zhh(1,1,kk);
        %             obj.ZdrProfileSim(kk,jj) = obj.ex.Zdr(1,1,kk);
        %         end
        %     end

        %     obj.ZhhProfileSim(imag(obj.ZhhProfileSim)>0) = -inf;
        %     obj.ZdrProfileSim(imag(obj.ZhhProfileSim)>0) = 0;
        % end

        %% use the forward model LUT to estimate the radar profiles from the DSD at the top of the rainshaft
        function varargout = estimateSingleRadarProfile(obj, N0,mu, lambda)
            % this is our backwards model that estimates the DSD profile.
            % in: Zhh in dBZ, Zdr in dB, Dm in mm for the top of the domain
            % out: dN, Zhhp, Zdrp, (dnint)

            %% calculate the DSD at the top of the rainshaft
            % keyboard
            Ntop = N0.*obj.ra.D.^(mu).*exp(-lambda.*obj.ra.D);
            Zdr = obj.ra.calcZdr(Ntop);
            Zhh = obj.ra.calcZhh(Ntop);
            [~,ind] = max(Ntop);
            Dm = obj.ra.D(ind);

            %% initialize the forward model interpolant
            % detaf = griddedInterpolant({obj.ZdrGrid, obj.RHGrid, obj.DmGrid}, obj.deta.(obj.bandName), 'linear', 'linear');

            if isempty(obj.detaInterpolant)
                if obj.hires
                    obj.detaInterpolant = load(sprintf('../data/LUTs/detaLUTs-N0MuLambda-hires-%1.0f.mat', obj.minuteIdx), 'detaInterpolant');
                else
                    obj.detaInterpolant = load(sprintf('../data/LUTs/detaLUTs-N0MuLambda-%1.0f.mat', obj.minuteIdx), 'detaInterpolant');
                end
                obj.detaInterpolant = obj.detaInterpolant.detaInterpolant;
            end

            %% calculate the difference kernel (detakk) at each height, RH
            detakk = squeeze(obj.detaInterpolant(...
                mu.*ones(size(obj.RHProfileObs)), ...
                lambda.*ones(size(obj.RHProfileObs)), ...
                obj.RHProfileObs ...
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

            Zhhp = (Zhhp)';
            Zdrp = (Zdrp)';

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

        %% initialize N based on reflectivity, Zdr, and Dm (not recommended)
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

        %% helper functions to convert between N0, mu, lambda and Zhh
        function N0 = getN0FromZhhMuLambda(obj, Zhh, mu, lambda)
            % in: Zhh in dBZ, mu, and lambda in mm^-1 
            % out: N0 in m^-3 mm^-1
            N1 = (obj.ra.D.^(mu) .* exp(-lambda.*obj.ra.D));
            Zhh1 = obj.ra.calcZhh(N1);
            N0 = 10.^(Zhh/10)./10.^(Zhh1/10);
        end

        %% helper function to use the gamma distribution to get N from N0, mu, lambda
        function [N] = getNFromN0MuLambda(obj, N0, mu, lambda)
            % in: N0 in m^-3 mm^-1, mu, and lambda mm^-1 
            % out: N in m^-3 mm^-1
            N = N0(:) .* (obj.ra.D.^(mu(:))) .* exp(-lambda(:).*obj.ra.D);
            N(isnan(N)) = 0;
        end

        %% estimate the N0, mu, and lambda from a given N using curve fitting.
        function [N0,mu,lambda] = getN0MuLambdaFromN(obj, N)
            % in: N in m^-3 mm^-1
            % out: N0 in m^-3 mm^-1, mu, and lambda in mm^-1 


            fun = @(x) calcErr(N, x(1), x(2), x(3));
            % initial guesses
            N01 = 1e8;
            mui = 5;
            lambdai = 5;
            % search for the optimal gamma and mu
            [x,~] =  fminsearchbnd(fun, [N01, mui, lambdai], [1, -2, 0], [1e11, 20, 20]);
            N0 = x(1);
            mu = x(2);
            lambda = x(3);

            

            function err = calcErr(N, N0, mu, lambda)
                Nsim = obj.getNFromN0MuLambda(N0, mu, lambda);
                % keyboard
                opt = optimset('TolX',1e-2,'TolFun',1e-2,'Display','off');
                err = rms((N(:) - Nsim(:))./max(N(:)), "omitmissing");
            end
        end
    end

    % methods 
    %     function bakeDSDs(obj)
    %         % function that calculates the DSD at the top of the rainshaft for many combinations of Zhh, Zdr, Dm
    %         % and saves them in a LUT. These can then be used later to interpolate the DSD at the top of the rainshaft for quick access.

    %         obj.ex = exmiras;

    %         NtopLUT = struct(...
    %             'S', NaN(numel(obj.ZhhGrid), numel(obj.ZdrGrid), numel(obj.DmGrid), 250), ...
    %             'C', NaN(numel(obj.ZhhGrid), numel(obj.ZdrGrid), numel(obj.DmGrid), 250), ...
    %             'X', NaN(numel(obj.ZhhGrid), numel(obj.ZdrGrid), numel(obj.DmGrid), 250) ...
    %         );
    %         progress = 0;
    %         total = numel(obj.ZhhGrid) * numel(obj.ZdrGrid) * numel(obj.DmGrid);
    %         for bandName = ["S", "C", "X"]             
    %             for Zhh = obj.ZhhGrid
    %                 for Dm = obj.DmGrid
    %                     for Zdr = obj.ZdrGrid
    %                         obj.ex = exmiras;

    %                         if mod(progress/total*100, 10) == 0
    %                             fprintf('Progress: %d/%d\n', progress, total)
    %                         end
            
    %                         % Run the ideal simulation                  
    %                         %%! set up initial droplet size distribution
    %                         % fprintf('T0m = %1.2f, humidity = %1.2f, dBZStart = %d, zdr = %1.2f\n', T0m, humidity, dBZStart, zdr)
    %                         obj.ex.rngToggle = false;
    %                         %% set up initial droplet size distribution
    %                         obj.ex.xgrid = [50];
    %                         obj.ex.ygrid = [50];
    %                         obj.ex.zgrid = 25:50:2025;

    %                         sx = numel(obj.ex.xgrid);
    %                         sy = numel(obj.ex.ygrid);
    %                         sz = numel(obj.ex.zgrid);
    %                         Zhhi = Zhh + zeros(sx, sy, sz)*0;
    %                         Zdri = Zdr + rand(sx, sy, sz)*0;
    %                         Zhhi(:,:,1:size(Zhhi,3)-1) = -inf;
    %                         Zdri(:,:,1:size(Zhhi,3)-1) = 0;

    %                         obj.ex = obj.ex.initFromLambdaName(bandName);
    %                         obj.ex = obj.ex.initFromDm(Zhhi(end), Zdri(end), Dm);
                            
    %                         NtopLUT.(bandName)(...
    %                             obj.ZhhGrid ==Zhh, ...
    %                             obj.ZdrGrid == Zdr, ...
    %                             obj.DmGrid == Dm, ...
    %                             : ...
    %                         ) = squeeze(obj.ex.N(1,1,end,:));

    %                         progress = progress + 1;
    %                     end
    %                 end
    %             end
    %         end

    %         ZhhGrid = obj.ZhhGrid;
    %         ZdrGrid = obj.ZdrGrid;
    %         DmGrid = obj.DmGrid;
    %         save('../data/LUTs/NtopLUTs.mat', 'NtopLUT', 'ZhhGrid', 'ZdrGrid', 'DmGrid')

    %     end
    % end

    %% data assimilation methods
    methods
        function varargout=profileOptimizer(obj, ZhhProfileObs, ZdrProfileObs, varargin)
            % function that optimizes a profile of simulated Zhh and Zdr to find the best fitting DSD.
            % inputs: ZhhProfileObs, ZdrProfileObs: observed profiles of Zhh and Zdr to fit
            % outputs: N0, mu, lambda: optimized values of N0, mu, and lambda
            % outputs: ~,  ~,  ~, fv: value of the objective function at the optimum
            % outputs: ~,  ~,  ~, ~, N: optimized DSD vertical profile

            options = optimset();
            options.TolX = 1e-1;
            options.TolFun = 0.5;

            %% parse optional inputs
            p = inputParser;
            addParameter(p, 'KdpProfileObs', [], @(x) isnumeric(x) || islogical(x));
            parse(p, varargin{:});
            KdpProfileObs = p.Results.KdpProfileObs;
        
            %% run the assimilation optimization routine
            % initial guess: median of the observed profiles, and 0.7 mm for D
            % use the Zhang et al. (2001) method to get an initial guess for N0, mu, and gamma, then run ensemble around that.
            lambdai = max(min(4./(mean(ZdrProfileObs(:), "omitmissing")),12), 0); % initial guess for gamma
            mui = -0.016*lambdai.^2 + 1.213*lambdai - 1.957;% initial guess for mu
            N0i = obj.getN0FromZhhMuLambda(mean(real(ZhhProfileObs(:)), "omitmissing"), mui, lambdai);
            x = [N0i, mui, lambdai];

            muList = mui + (-2:1:2);
            lambdaList = lambdai + (-2:1:2);
            muList = muList(muList>=-2 & muList<=15);
            lambdaList = lambdaList(lambdaList>=0 & lambdaList<=20);
            try
                for ii= 1:length(muList)
                    for jj = 1:length(lambdaList)
                        xStart = [N0i, muList(ii), lambdaList(jj)];
                        [x,fv]=fminsearchbnd(@errorFunc, ...
                            xStart, ...
                            [1, -2, 0], ...
                            [1e11, 15, 20], ...
                            options);
                        xtemp(ii,:) = x;
                        fvtemp(ii) = fv;
                        [fv, ind] = min(fvtemp);
                        x = xtemp(ind, :);
                    end     
                end
            catch
                warning('Optimization failed, returning NaNs')
                switch nargout
                    case 3
                        varargout = {NaN, NaN, NaN};
                    case 4
                        varargout = {NaN, NaN, NaN, NaN};
                    case 5
                        varargout = {NaN, NaN, NaN, NaN, NaN};
                end
                return
            end 

            %% extract the optimized and return outputs.
            N0Opt = x(1); % dBZ
            muOpt = x(2); % dB
            lambdaOpt = x(3); % mm
            [N,~,~] = obj.estimateSingleRadarProfile(N0Opt, muOpt, lambdaOpt);
            switch nargout
                case 3
                    varargout = {N0Opt, muOpt, lambdaOpt};
                case 4
                    varargout = {N0Opt, muOpt, lambdaOpt, fv};
                case 5
                    % keyboard
                    varargout = {N0Opt, muOpt, lambdaOpt, fv, N};
            end

            function errorReturn=errorFunc(x)
                % calculates the error between the observed and simulated profiles for a given set of Zhh, Zdr, and Dm
                % x(1) = N0, x(2) = mu, x(3) = lambda
                if isa(x, 'table')
                    x = table2array(x);
                end
            
                N0 = x(1); % mm^-1 m^-3
                mu = x(2); % 
                lambda = x(3); % mm^-1
                
                
                % function to optimize the profile of the DSD
                [~, Zhhp, Zdrp,Kdpp] = obj.estimateSingleRadarProfile(N0, mu, lambda);
                
                % solve for Zvv from Zhh and Zdr
                Zvvp = Zhhp - Zdrp;
                
                % convert to linear units
                ZhhProfileObsLinear = 10.^(ZhhProfileObs./10);
                Zhhp = 10.^(Zhhp./10);
                Zvvp = 10.^(Zvvp./10);
                
                %% calculate the error between the observed and simulated profiles
                % Zhh
                errorRawZhh = (Zhhp(:) - ZhhProfileObsLinear);
                errorZhh = rms(mean(errorRawZhh, 'all', 'omitnan'), "omitmissing");%./std(ZhhProfileObsLinear(:), "omitmissing");
                
                %Zdr
                errorRawZvv = (Zvvp(:) - 10.^((ZhhProfileObs - ZdrProfileObs)./10));
                errorZvv = rms(mean(errorRawZvv, 'all', 'omitnan'), "omitmissing");%./std(Zvvp(:), "omitmissing");

                % Kdp
                
                if ~isempty(KdpProfileObs)
                    % keyboard
                    KdpProfileObsLinear = (KdpProfileObs);
                    errorRawKdp = (Kdpp(:) - KdpProfileObsLinear);
                    errorKdp = rms(mean(errorRawKdp, 'all', 'omitnan'), "omitmissing")./0.05;
                else
                    errorKdp = 0;
                end

                errorReturn = (mean([errorZhh, errorZvv, errorKdp], "all","omitnan")); % + errorModeZdr.^2;

                % penalize solutions that are outside the physical bounds of the forward model
                polytop = [2.3277, -0.9189];
                polybot = [0.2796, -2.0];
                distbot = x(2) - polyval(polybot, x(3));
                disttop = x(2) - polyval(polytop, x(3));
                if distbot < 0 
                    errorReturn = errorReturn + (distbot).^2;
                end

                if disttop > 0
                    errorReturn = errorReturn + (disttop).^2;
                end

                if isnan(errorReturn)
                    % keyboard
                    error('not enough data to constrain the optimization');
                end
                

            end
        end

        %% alias of pointOptimizer. kept to remain consistent with profileOptimizer()
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

    %% constructor
    methods 
        function da = dsdAssimilation(lambda)
            % loads an instance of the radar class for the specified wavelength
            da.ra = radar(lambda);
            da.De = logspace(log10(0.1),log10(da.dMax),da.nBins+1);
            da.D =  (da.De(1:end-1) + da.De(2:end))/2;
            da.Dw = diff(da.De);

        end
    end
end
