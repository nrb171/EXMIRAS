classdef exmiras < thermo
    % evaporation model, hooked into thermo class.
    %% set properties
    properties
        
        pv % hPa; partial pressure of water vapor

        %% kinematic properties
        u % zonal (0,1,0) wind speed in m/s (same size as T)
        v % meridional (1,0,0) wind speed in m/s (same size as T)
        w % vertical (0,0,1) wind speed in m/s (same size as T)


        %% thermodynamic properties
        dHevap  = 2256e3 % J/kg; enthalpy of vaporization


        N       % m^-3 mm^-1; number concentration of droplets
        N0      % m^-3; previous number concentration of droplets
        N00    % initial scaling factor for number concentration of droplets
        NP      % m^-3; number for "parent" field (nonzero where cloud droplets are constantly being replenished and reset at each loop of integration)
        mu      %shape parameter
        gamma   % shape parameter

        nBins = 250 %number of bins for droplet size distribution
        dt = 10 % s; time step
        st % length of storm (raining) in seconds

        %% simulation domain (same size as T)
        xgrid
        ygrid
        zgrid
        
        
        lambda = 111% mm/band name; wavelength of radar
        dpp % dual-pol preprocessing structure

        De % mm; edges of droplet size bins
        Dw % mm; width of droplet size bins
        D % mm; characteristic droplet size in each bin
        M % kg; mass of droplets in each bin
        vt % function of D; terminal velocity of droplets
        



    end

    properties (Access = private)
        nSteps = 0 % number of steps taken in integration

        epsilon = 62.1 - sqrt(-1)*32.0; % dielectric constant of liquid water
        
        Kw = 0.93; % reflectivity factor of water

        refCalc = 'NEW'; % method of calculating reflectivity

        E % collision efficiency
        
    end

    %% constant properties
    properties (Constant, Access = private)
        dMax = 8 % mm; %maximum droplet size 

    end

    properties (Dependent)

        m % kg/m^3 total mass of liquid water in air volume
        mv % kg/m^3 total mass of water vapor in air volume
        theta % potential temperature


      

        %% diagnostic quantities
        Zhh % reflectivity (horizontal)
        Zvv % reflectivity (vertical)
        Zdr % differential reflectivity
        rhohv % correlation coefficient
        kdp % specific differential phase shift
        RR % rain rate

        %% rates of change
        dmevap % change in mass due to evaporation
        dDevap % change in droplet size due to evaporation
        dNevap % change in number concentration due to evaporation

        dNadvect % change in number concentration due to advection
        dNfallout % change in number concentration due to fallout
        dNcoal % change in number concentration due to collision-coalescence
        
        dTevap % change in temperature from evaporative process

        vol % m^3; volume of each grid cell
        qv %
        
    end

    
    methods
        function obj = initFromLambdaName(obj, lambdaName)
            % warning('off')
            
                bandName = ["S",    "C",    "X",        "Ku",       "Ka"       ];
                wavelength = [111,  53.5,   33.3,       22,         8.43       ];
                obj.lambda = wavelength(strcmp(bandName, lambdaName));
            % end
            % warning('on')
        end
        
        
        % change in droplet size due to evaporation
        function dD = get.dDevap(obj)
            % set up change in droplet size
            S = obj.pv./obj.es;
            r = obj.D/1000/2;
            gamma = ((S - 1)./(obj.FK + obj.FD));
            gamma = reshape(gamma, [], obj.nBins);
            drdt = gamma./r;
            
            dr = drdt*obj.dt;
            dD = dr*2*1000;
            dD = reshape(dD, [size(obj.T), obj.nBins]);
        end

        function mv = get.mv(obj)
            % out: mass of water vapor per unit volume in kg/m^3
            mv = 100*obj.pv./(obj.Rv.*obj.T);
        end

        function theta = get.theta(obj)
            theta = obj.T.*(1000./obj.p).^(obj.Ra/obj.cp);
        end

        function qv = get.qv(obj)
            qv = obj.pv./obj.es;
        end
        

    end

    %% reflectivity methods
    methods
        function dpp = dualPolPreprocessing(obj, lambda)
            % various preprocessing steps to calculate dual-pol variables
            % keyboard
            r = 0.9951 + 0.02510*obj.D - 0.03644*obj.D.^2 + 0.005303*obj.D.^3 - 0.0002492*obj.D.^4;
            b = ((obj.D/1000).^3./r).^(1/3);
            a = b.*r;
            % keyboard
            
            S = load('S.mat');
            iWavelength = find(S.wavelengths == lambda);
            % keyboard

            dpp.fa = cellfun(@(x) x(1,1), S.S(iWavelength,:));
            dpp.fb = -cellfun(@(x) x(2,2), S.S(iWavelength,:));

            reflh = readmatrix('./LUTs/reflh.txt');
            reflv = readmatrix('./LUTs/reflv.txt');
            xsecth = readmatrix('./LUTs/xsecth.txt');
            xsectv = readmatrix('./LUTs/xsectv.txt');
            rhohva = readmatrix('./LUTs/rhohva.txt');
            rhohvb = readmatrix('./LUTs/rhohvb.txt');
            kdpl = readmatrix('./LUTs/kdp.txt');
            % N = squeeze(ev.N(1,1,end,:))';


            iWavelength = find([111, 53.5, 33.3, 22.0, 8.43] == lambda);
            dpp.reflh = reflh(iWavelength,:);
            dpp.reflv = reflv(iWavelength,:);
            dpp.xsecth = xsecth(iWavelength,:);
            dpp.xsectv = xsectv(iWavelength,:);
            dpp.rhohva = rhohva(iWavelength,:);
            dpp.rhohvb = rhohvb(iWavelength,:);
            dpp.kdpl = kdpl(iWavelength,:);
            
            r1 = exp(-2*(10*1/180*pi)^2);
            dpp.A1 = 1/4*(1+r1)^2;
            dpp.A2 = 1/4*(1-r1)^2;
            dpp.A3 = (3/8 + 1/2*r1+1/8*r1^4)^2;
            dpp.A4 = (3/8 - 1/2*r1+1/8*r1^4)*(3/8 + 1/2*r1+1/8*r1^4);
        end
        %% calculators for reflectivity
        function Zhh = calcZhh(obj, N, dpp)
            % in: N in m^-3 mm^-1, D in mm
            % out: horizontal reflectivity in dBZ
            % keyboard
            N = reshape(N, [], obj.nBins);
        

            
            % dpp = obj.dualPolPreprocessing(obj.lambda);
            if strcmp(obj.refCalc, 'OLD')
                Zhh = 10*log10(...
                    4*obj.lambda.^4./(pi^4*abs(obj.Kw).^2)*...
                    trapz(obj.D,...
                        (abs(dpp.fb).^2 - 2*real(conj(dpp.fb).*(dpp.fb-dpp.fa)).*dpp.A2 + abs(dpp.fb-dpp.fa).^2.*dpp.A4).*N, 2)...
                );
            else
                Zhh = 10*log10(...
                    trapz(obj.D, dpp.reflh.*N, 2) ...
                );
            end
            if numel(N) == obj.nBins
                Zhh = Zhh;
            else
                Zhh = reshape(Zhh, numel(obj.xgrid), numel(obj.ygrid), numel(obj.zgrid));
            end
        end

        function Zvv = calcZvv(obj, N, dpp)
            % in: N in m^-3 mm^-1, D in mm
            % out: vertical reflectivity in dBZ
            N = reshape(N, [], obj.nBins);
            % dpp = obj.dualPolPreprocessing(obj.lambda);
            if strcmp(obj.refCalc, 'OLD')
                Zvv = 10*log10(4*obj.lambda.^4./(pi^4*abs(obj.Kw).^2)*trapz(obj.D,(abs(dpp.fb).^2 - 2*real(conj(dpp.fb).*(dpp.fb-dpp.fa)).*dpp.A1 + abs(dpp.fb-dpp.fa).^2.*dpp.A3).*N, 2));    
            else
                Zvv = 10*log10(...
                    trapz(obj.D, dpp.reflv.*N, 2) ...
                );
            end
            if numel(N) == obj.nBins
                Zvv = Zvv;
            else
                Zvv = reshape(Zvv, numel(obj.xgrid), numel(obj.ygrid), numel(obj.zgrid));
            end
        end

        function rhohv = calcRhohv(obj, N, dpp)
            N = reshape(N, [], obj.nBins);
            % D = obj.D;
            % iWavelength = 1
            a = dpp.rhohva;
            b = dpp.rhohvb;
            c = dpp.xsecth./(2*pi);
            d = dpp.xsectv./(2*pi);
            % keyboard
            % rhohv = readmatrix('/h/eol/nbarron/workshop/apar-scripts/evapModel/LUTs/rhohv.txt');
            % N([1,3:end]) = 0;
            % rho_hv = trapz(D, N.*rhohv(iWavelength,:))/trapz(D,N);
            % N = N(end,:)
            rho_hv = ... 
                sqrt(...
                    ((trapz(obj.D, N.*a,2)).^2 + ...
                    (trapz(obj.D, N.*b,2)).^2)) ./ ...
                    sqrt((trapz(obj.D, N.*c,2) .* trapz(obj.D, N.*d,2)) ...
            );
            rhohv = reshape(rho_hv, size(obj.N, [1:3]));
        end
            
        function kdp = calcKdp(obj, N, dpp)
            N = reshape(N, [], obj.nBins);
            kdp = trapz(obj.D, N.*dpp.kdpl,2);
            kdp = reshape(kdp, size(obj.N, [1:3]));
        end
        %% initialize N based on reflectivity and Zdr
        function obj = initFromReflectivity(obj, Z, Zdr)
            % in: Z in dBZ, Zdr in dBz
            % out: N in m^-3 mm^-1
            % keyboard
            inds = find(~isinf(Z));
            [ix, iy, iz] = ind2sub(size(Z), inds);
            mu = NaN(size(Z));
            gamma = NaN(size(Z));
            % keyboard
            obj.dpp = obj.dualPolPreprocessing(obj.lambda);
            for i = 1:length(inds)
                ind = inds(i);
                % keyboard
                fun = @(x) calcErr(obj, obj.dpp, Z(ind), Zdr(ind), initN(obj, obj.dpp, Z(ind), x(1), x(2)), x(1), x(2));
                [x,fv] =  fminsearchbnd(fun, [5, 2.5], [0, -1], [20, 5]);

                obj.N(ix(i), iy(i), iz(i), :) = initN(obj, obj.dpp, Z(ind), x(1), x(2));
                mu(ix(i), iy(i), iz(i)) = x(2);
                gamma(ix(i), iy(i), iz(i)) = x(1);
            end
            % keyboard
            obj.mu = mu;
            obj.gamma = gamma;
            obj.N00 = getN0(obj, obj.dpp, Z, gamma, mu);

            function N0 = getN0(obj, dpp, dBZi, gamma, mu)
                % various preprocessing steps to calculate dual-pol variables
                
                % keyboard
                
                mu = mu(:);
                gamma = gamma(:);


                %! overriding mu
                % mu = zeros(size(gamma));
                % mu(isnan(gamma)) = NaN;
                % mu = -0.0201*gamma.^2 + 0.902*gamma - 1.78;


                % N = reshape(N, [], obj.nBins);
                N1 = (obj.D.^(mu) .* exp(-gamma.*obj.D));
                % N1 = repmat(N1, [1, 1])
                Zhh = calcZhh(obj, N1, dpp);

                N0 = 10.^(dBZi(:)/10)./10.^(Zhh(:)/10);
                N0 = reshape(N0, size(obj.mu));
            end

            function [N] = initN(obj, dpp, dBZi, gamma, mu)
                % various preprocessing steps to calculate dual-pol variables

                % N = reshape(N, [], obj.nBins);
                N1 = (obj.D.^(mu) .* exp(-gamma*obj.D));
                % N1 = repmat(N1, [1, 1])
                Zhh = 10^(calcZhh(obj, N1, dpp)/10);
                % Zvv = calcZvv(obj, N1, dpp);
                % Zhh = 4*obj.lambda.^4./(pi^4*abs(obj.Kw).^2)*trapz(obj.D,(abs(dpp.fb).^2 - 2*real(conj(dpp.fb).*(dpp.fb-dpp.fa)).*dpp.A2 + abs(dpp.fb-dpp.fa).^2.*dpp.A4).*N1);
                % Zvv = 4*obj.lambda.^4./(pi^4*abs(obj.Kw).^2)*trapz(obj.D,(abs(dpp.fb).^2 - 2*real(conj(dpp.fb).*(dpp.fb-dpp.fa)).*dpp.A1 + abs(dpp.fb-dpp.fa).^2.*dpp.A3).*N1);
                % keyboard

                

                N0 = 10^(dBZi/10)./(Zhh);

                N = N0(:) .* (obj.D.^(mu) .* exp(-gamma*obj.D));
            end
            function err = calcErr(obj, dpp, dBZi, Zdri, N, gamma, mu)
                % err = 1;
                % keyboard
                % b = obj.D;
                % r = 0.9951 + 0.02510*obj.D - 0.03644*obj.D.^2 + 0.005303*obj.D.^3 - 0.0002492*obj.D.^4;
                % a = obj.D.*r;

                % obj.D = (a.*b.^2).^(1/3);
                Zhh2 = calcZhh(obj, N, dpp);
                Zvv2 = calcZvv(obj, N, dpp);
                % Zhh2 = 10*log10(4*obj.lambda.^4./(pi^4*abs(obj.Kw).^2)*trapz(obj.D,(abs(dpp.fb).^2 ).*N, 2))
                % Zvv2 = 10*log10(4*obj.lambda.^4./(pi^4*abs(obj.Kw).^2)*trapz(obj.D,(abs(dpp.fa).^2 ).*N, 2))
                % Zdr = Zhh2 - Zvv2;

                % keyboard
                dZdr = Zhh2 - Zvv2 - Zdri;
                err = abs(dZdr);
                
            end
        end

        %% get methods for reflectivity
        function Zhh = get.Zhh(obj)
            % in: D in mm
            % out: horizontal reflectivity in dBZ
            % keyboard
            % dpp = obj.dualPolPreprocessing(obj.lambda);
            Zhh = obj.calcZhh(obj.N, obj.dpp);
        end
        function Zvv = get.Zvv(obj)
            % in: D in mm
            % out: vertical reflectivity in dBZ
            % dpp = obj.dualPolPreprocessing(obj.lambda);
            Zvv = obj.calcZvv(obj.N, obj.dpp);
        end
        function Zdr = get.Zdr(obj)
            Zdr = obj.Zhh - obj.Zvv;
        end
        function rhohv = get.rhohv(obj)
            % dpp = obj.dualPolPreprocessing(obj.lambda);
            rhohv = obj.calcRhohv(obj.N, obj.dpp);
        end
        function kdp = get.kdp(obj)
            % dpp = obj.dualPolPreprocessing(obj.lambda);
            kdp = obj.calcKdp(obj.N, obj.dpp);
        end

        %% rain rate based on reflectivity
        function RR = get.RR(obj)
            % in: D in mm
            % out: RR in mm/hr
            N = reshape(obj.N, [], obj.nBins);

            % dpp = obj.dualPolPreprocessing(obj.lambda);
            RR = pi/6*trapz(obj.D, (obj.D/1000).^3 .* N .* obj.vt, 2)*3600*1000;
            RR = reshape(RR, size(obj.N, [1:3]));
        end 
    end
    %% Model core
    methods 
        % change in dsd due to evaporation
        function dN = get.dNevap(obj)
            % set up change in number concentration            
            dNt = reshape(obj.N, [], obj.nBins).*(reshape(obj.dDevap, [], obj.nBins) ./ (obj.Dw));            
            dN = dNt;
            
            dN(:,1:end-1) = dN(:,1:end-1)-dNt(:,2:end);
            dN = reshape(dN, [size(obj.T), obj.nBins]);
        end

        % function dN = get.dNadvect(obj)
            
        % end

        function dN = get.dNfallout(obj)
            %% get change in number concentration due to fallout
            % calculate the 'true' vertical motion of droplets by subtracting terminal velocity from vertical velocity
            % keyboard
            wTrue = reshape(obj.w(:) - obj.vt, [size(obj.w), length(obj.vt)]);
            dz = wTrue*obj.dt;
            
            % keyboard

            %% calculate fraction of droplets that fall out/pushed upward
            % positive values indicate upward motion, negative values indicate downward motion
            r = dz./repmat(permute(repmat(gradient(obj.zgrid), [numel(obj.xgrid),1,numel(obj.ygrid)]), [3,1,2]), [ones(1,numel(size(obj.w))),length(obj.vt)]);

            %% find the indices of upward/downward motion, then add/subtract one from the appropriate indices.
            % upward motion
            indsUp0 = find(r > 0);
            [iUp,jUp,kUp,lUp] = ind2sub(size(r), indsUp0);
            kUp = kUp + 1;
            maskUp = kUp > size(obj.N,3);
            iUp(maskUp) = [];
            jUp(maskUp) = [];
            kUp(maskUp) = [];
            lUp(maskUp) = [];
            indsUp0(maskUp) = [];
            indsUp = sub2ind(size(r), iUp, jUp, kUp, lUp);

            % downward motion
            indsDown0 = find(r < 0);
            [iDown,jDown,kDown,lDown] = ind2sub(size(r), indsDown0);
            kDown = kDown -1;
            maskDown = kDown < 1;
            iDown(maskDown) = [];
            jDown(maskDown) = [];
            kDown(maskDown) = [];
            lDown(maskDown) = [];
            indsDown0(maskDown) = [];
            indsDown = sub2ind(size(obj.N), iDown, jDown, kDown, lDown);


            %% calculate change in number concentration due to fallout/positive vertical advection
            dN = -abs(obj.N.*r);

            dN(indsUp) = dN(indsUp) + abs(obj.N(indsUp0).*r(indsUp0));
            dN(indsDown) = dN(indsDown) + abs(obj.N(indsDown0).*r(indsDown0));

            
        end

        function E=getE(obj)
            % function to initialize collision efficiency LUT (used in constructor)

            vi = repmat(obj.vt', [1,size(obj.N,4)]);
            vj = repmat(obj.vt, [size(obj.N,4),1]);
            di = repmat(obj.D', [1,size(obj.N,4)])/1000;
            dj = repmat(obj.D, [size(obj.N,4),1])/1000;

            sigma = 7.28e-2;
            a = 0.778;
            b = 2.61e6;
            CKE = obj.rhol*pi/12.*dj.^3.*di.^3./(di.^3 + dj.^3).*(vj-vi).^2;
            SC = pi * sigma * (dj.^3 + di.^3).^(2/3);
            ST = pi * sigma * (dj.^2 + di.^2);
            deltaS = ST - SC;
            ET = CKE + deltaS;

            E = a*(1+di./dj).^(-2).*exp(-b*sigma * ET.^2./SC);
            E(ET>5e-6) = 0;
        end
        function dN = get.dNcoal(obj)

            dN = zeros(size(obj.N));
            inds = find(sum(obj.N, 4)~=0);
            for ii = 1:numel(inds)
                ind = inds(ii);
                [i,j,k] = ind2sub(size(obj.N, [1:3]), ind);
                % keyboard
                % mi = repmat(obj.M', [1,size(obj.N,4)]);
                % species j collects species i. Integrate along columns
                % i.e. along the rows are the donor droplets, along the columns are the absorbing droplets
                vi = repmat(obj.vt', [1,size(obj.N,4)]);
                vj = repmat(obj.vt, [size(obj.N,4),1]);
                di = repmat(obj.D', [1,size(obj.N,4)])/1000;
                dj = repmat(obj.D, [size(obj.N,4),1])/1000;
                Ni = repmat(squeeze(obj.N(i,j,k,:)), [1,size(obj.N,4)]);
                Nj = repmat(squeeze(obj.N(i,j,k,:))', [size(obj.N,4),1]);

                % K = pi*E(1:149,150)'.*(obj.D(150)/2 + obj.D(1:149)/2).^2.*abs(obj.vt(150)-obj.vt(1:149));
                dmij = triu(obj.rhol*4*pi/3* ... %kg/
                    pi*((di+dj)/2).^2.*abs(vj-vi).* ... /s
                    (di/2).^3.*Ni.*Nj.*obj.E);

                % keyboard
                dmj = trapz(obj.D, dmij, 1); % mass gained per droplet 
                dmi = trapz(obj.D, dmij, 2); %mass lost per droplet bin 
                dmi = dmi * sum(dmj)./sum(dmi); % mass lost per droplet bin

                % keyboard
                % fig = figure("Units", "inches", "Position", [0,0,3.75,3]);
                % plot(obj.D, squeeze(dmj)./obj.M*obj.dt, 'r')
                % hold on 
                % plot(obj.D, -squeeze(dmi')./obj.M*obj.dt, 'b')
                % legend('number concentration gained', 'number concentration lost', 'Location', 'southeast')
                % xlabel('D [mm]')
                % ylabel('dN [m^{-3} mm^{-1}]')
                % xscale('log')
                % print2(fig, './dNcoal.png')
                dN(i,j,k,:) = (dmj-dmi')./obj.M*obj.dt; % number concentration gained per droplet bin
            end

            
        end


        function vol = get.vol(obj)
            dy = repmat(gradient(obj.ygrid), [numel(obj.xgrid), 1, numel(obj.zgrid)]);

            dx = repmat(gradient(obj.xgrid)', [1, numel(obj.ygrid), numel(obj.zgrid)]);
            dz = repmat(gradient(obj.zgrid), [numel(obj.xgrid), 1,numel(obj.ygrid)]);
            dz = permute(dz, [1,3,2]);
            vol = dx.*dy.*dz;
        end

        function dm = get.dmevap(obj)
            
            dNevapf = reshape(obj.dNevap, [], obj.nBins);
            
            dm = dNevapf .* obj.M;
            dm = reshape(dm, [size(obj.N)]);
            dm = trapz(obj.D,dm, numel(size(dm)));
            % dDevap = reshape(obj.dDevap, [], obj.nBins);
            % D0 = repmat(obj.D, [size(dDevap,1),1]);
            % D1 = D0 + dDevap;

            % v0 = 4/3*pi*(D0/2/1000).^3; 
            % v1 = 4/3*pi*(D1/2/1000).^3;
            % m0 = obj.rhol*v0;
            % m1 = obj.rhol*v1;

            % dm = obj.dNevap.*reshape((m0 - m1), size(obj.dDevap));
            % dm = sum(dm, numel(size(dm)));

        end

        function dT = get.dTevap(obj)    
            dQ = obj.dHevap*obj.dmevap;
            dT = dQ./(obj.rhoa*obj.cp); 
        end

        
    end

    %% constructor
    methods
        function obj = integrate(obj)
            if obj.nSteps == 0
                % keyboard
                if isempty(obj.w)
                    obj.w = zeros(size(obj.T));
                end
                if isempty(obj.u)
                    obj.u = zeros(size(obj.T));
                end
                if isempty(obj.v)
                    obj.v = zeros(size(obj.T));
                end
            end

            %% calculate rates of change
            % keyboard
            dT = obj.dTevap;
            dN = obj.dNevap + obj.dNfallout + obj.dNcoal;
            dm = obj.dmevap;

            if any(obj.T(:)<abs(dT(:))) || any(obj.mv(:) < abs(dm(:)))
                obj.dt = obj.dt/2;
                warning('Time step too large, reducing to %d', obj.dt)
                return

            end

            %% update state variables

            obj.N0 = obj.N;
            obj.N = obj.N + dN;
            obj.T = obj.T + dT;
            obj.pv = (obj.mv - dm).*obj.Rv.*obj.T/100;

            %% update the parent field
            if obj.nSteps < obj.st
                obj.N(obj.NP > 0) = obj.NP(obj.NP > 0);
            end

            %% update the number of steps
            obj.nSteps = obj.nSteps + 1;
        end

        function obj = evap(obj)

            % calculate dropsize distribution
            obj.De = logspace(log10(0.1),log10(obj.dMax),obj.nBins+1);
            obj.D =  (obj.De(1:end-1) + obj.De(2:end))/2;
            obj.Dw = diff(obj.De);

            alpha = 3.78; %m/s/mm^0.67
            beta = 0.67;
            %vt = alpha*obj.D.^beta(rho0/obj.rhoa)^0.4;
            obj.vt = -0.1021 + 4.932*obj.D - 0.955*obj.D.^2 + 0.07934*obj.D.^3 - 0.002362*obj.D.^4;

            obj.M = (4*pi/3*obj.rhol*(obj.D/1000/2).^3);

            obj.E = obj.getE();
            
        end

    end
    
end