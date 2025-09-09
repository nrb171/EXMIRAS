classdef exmiras < thermo
    % evaporation model, hooked into thermo class.
    %% set properties
    properties
        
        pv % hPa; partial pressure of water vapor

        %% kinematic properties
        u % zonal (0,1,0) wind speed in m/s (same size as T)
        v % meridional (1,0,0) wind speed in m/s (same size as T)
        w % vertical (0,0,1) wind speed in m/s (same size as T)


        

        %% number concentration properties
        N       % m^-3 mm^-1; number concentration of droplets
        N0      % m^-3; previous number concentration of droplets
        N00    % initial scaling factor for number concentration of droplets
        NP      % m^-3; number for "parent" field (nonzero where cloud droplets are constantly being replenished and reset at each loop of integration)
        mu      %shape parameter
        gamma   % shape parameter

        dt = 10 % s; time step
        st % length of simulation in seconds

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

        ts % time coordinate for the simulation, in seconds

        waterVaporInteraction = true % whether to include water vapor interaction (evaporation) or not
        coalToggle = true % whether to include coalescence or not
        falloutToggle = true % whether to include fallout or not
        evapToggle = true % whether to include evaporation or not

        rngToggle = true

        %% assistant classes
        ra  % radar class instance
        da % dsdAssimilation class instance

        



    end

    properties (Access = private)
        nSteps = 0 % number of steps taken in integration

        epsilon = 62.1 - sqrt(-1)*32.0; % dielectric constant of liquid water
        
        Kw = 0.93; % reflectivity factor of water

        % refCalc = 'NEW'; % method of calculating reflectivity

        E % collision efficiency
        
    end

    %% constant properties
    properties (Constant)
        nBins = 250 %number of bins for droplet size distribution

        dMax = 8 % mm; %maximum droplet size 

        %% thermodynamic properties
        dHevap  = 2256e3 % J/kg; enthalpy of vaporization
    end

    properties (Dependent)

        m % kg/m^3 total mass of liquid water in air volume
        mv % kg/m^3 total mass of water vapor in air volume
        theta % potential temperature


      

        %% radar quantities
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
            dD = min(dD, repmat(obj.D, [size(dD,1), 1])); % limit to maximum droplet size
            % dD = 
            
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

    %% radar methods
    methods
        
        function obj = initializeRadarSimulator(obj, lambdaName)
            %% initialize the radar object from a band name
            bandName = ["S",    "C",    "X",        "Ku",       "Ka"       ];
            wavelength = [111,  53.5,   33.3,       22,         8.43       ];
            obj.lambda = wavelength(strcmp(bandName, lambdaName));

            obj.ra = radar(obj.lambda);
            obj.ra.rngToggle = obj.rngToggle;

            obj.da = dsdAssimilation(obj.lambda);
            obj.da.ra.rngToggle = obj.rngToggle;

        end

        %% get methods for reflectivity
        function Zhh = get.Zhh(obj)
            % in: D in mm
            % out: horizontal reflectivity in dBZ
            % keyboard
            % dpp = obj.dualPolPreprocessing(obj.lambda);
            Zhh = zeros(size(obj.N, [1:3]));
            inds = find(nansum(obj.N,4));
            for i = 1:length(inds)
                ind = inds(i);
                [xi, yi, zi] = ind2sub(size(obj.N, [1:3]), ind);
                
                Zhh(ind) = obj.ra.calcZhh(squeeze(obj.N(xi,yi,zi,:)));
            end
            % Zhh = obj.ra.calcZhh(obj.N);


        end
        function Zvv = get.Zvv(obj)
            % in: D in mm
            % out: vertical reflectivity in dBZ
            Zvv = zeros(size(obj.N, [1:3]));
            inds = find(nansum(obj.N,4));
            for i = 1:length(inds)
                ind = inds(i);
                [xi, yi, zi] = ind2sub(size(obj.N, [1:3]), ind);
                Zvv(ind) = obj.ra.calcZvv(squeeze(obj.N(xi,yi,zi,:)));
            end
        end
        function Zdr = get.Zdr(obj)
            Zdr = obj.Zhh - obj.Zvv;
        end
        function rhohv = get.rhohv(obj)
            rhohv = zeros(size(obj.N, [1:3]));
            inds = find(nansum(obj.N,4));
            for i = 1:length(inds)
                ind = inds(i);
                [xi, yi, zi] = ind2sub(size(obj.N, [1:3]), ind);
                rhohv(ind) = obj.ra.calcRhohv(squeeze(obj.N(xi,yi,zi,:)));
            end
        end
        function kdp = get.kdp(obj)
            kdp = zeros(size(obj.N, [1:3]));
            inds = find(nansum(obj.N,4));
            for i = 1:length(inds)
                ind = inds(i);
                [xi, yi, zi] = ind2sub(size(obj.N, [1:3]), ind);
                kdp(ind) = obj.ra.calcKdp(squeeze(obj.N(xi,yi,zi,:)));
            end
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
         
            dN(:,1:end-1) = (dN(:,1:end-1)-dNt(:,2:end));
            dN = reshape(dN, [size(obj.T), obj.nBins]);

            dN(obj.N+dN<0) = -obj.N(obj.N+dN<0); % fix overflowing values.
         
            inds = find((nansum(obj.N,4)));
            [xi, yi, zi] = ind2sub(size(obj.T), inds);
         
            % for ii = 1:numel(inds)
            %     dN(xi(ii), yi(ii), zi(ii),:) = filloutliers(squeeze(dN(xi(ii),yi(ii),zi(ii),:)),'linear','movmedian', 3);
            % end
         
        end

        % function dN = get.dNadvect(obj)
            
        % end

        function dN = get.dNfallout(obj)
            %% get change in number concentration due to fallout
            % calculate the 'true' vertical motion of droplets by subtracting terminal velocity from vertical velocity
            wTrue = reshape(obj.w(:) - obj.vt, [size(obj.w), length(obj.vt)]);
            dz = wTrue*obj.dt;
            
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

                % species j collects species i. Integrate along columns
                % i.e. along the rows are the donor droplets, along the columns are the absorbing droplets
                vi = repmat(obj.vt', [1,size(obj.N,4)]);
                vj = repmat(obj.vt, [size(obj.N,4),1]);
                di = repmat(obj.D', [1,size(obj.N,4)])/1000;
                dj = repmat(obj.D, [size(obj.N,4),1])/1000;
                Ni = repmat(squeeze(obj.N(i,j,k,:)), [1,size(obj.N,4)]);
                Nj = repmat(squeeze(obj.N(i,j,k,:))', [size(obj.N,4),1]);

                dmij = triu(obj.rhol*4*pi/3* ... %kg/
                    pi*((di+dj)/2).^2.*abs(vj-vi).* ... /s
                    (di/2).^3.*Ni.*Nj.*obj.E);
                dmj = trapz(obj.D, dmij, 1); % mass gained per droplet 
                dmi = trapz(obj.D, dmij, 2); %mass lost per droplet bin 
                dmi = dmi * sum(dmj)./sum(dmi); % mass lost per droplet bin

                dN(i,j,k,:) = (dmj-dmi')./obj.M*obj.dt; % number concentration gained per droplet bin
            end

            
        end


        function vol = get.vol(obj)
            % get the volume of each grid cell
            dy = repmat(gradient(obj.ygrid), [numel(obj.xgrid), 1, numel(obj.zgrid)]);

            dx = repmat(gradient(obj.xgrid)', [1, numel(obj.ygrid), numel(obj.zgrid)]);
            dz = repmat(gradient(obj.zgrid), [numel(obj.xgrid), 1,numel(obj.ygrid)]);
            dz = permute(dz, [1,3,2]);
            vol = dx.*dy.*dz;
        end

        function dm = get.dmevap(obj)
            % get change in mass due to evaporation
            
            dNevapf = reshape(obj.dNevap, [], obj.nBins);
            
            dm = dNevapf .* obj.M;
            dm = reshape(dm, [size(obj.N)]);
            dm = trapz(obj.D,dm, numel(size(dm)));
        end

        function dT = get.dTevap(obj)    
            dQ = obj.dHevap*obj.dmevap;
            dT = dQ./(obj.rhoa*obj.cp); 
        end
    end

    %% constructor
    methods
        function obj = integrate(obj)
            %% integrate one time step of the model
            % assumes that the state variables (T, p, pv, N, u, v, w) are already initialized
            % assumes that the dsd (N) is already initialized
            % iteratively calls the dN methods and updates the state variables
            

            % initialize wind fields if not provided
            if obj.nSteps == 0
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
            operators = ["dNevap", "dNfallout", "dNcoal"];
            
            % calculate dT and dm first, saving values for later 
            dT = obj.dTevap;
            dm = obj.dmevap;

            % record previous N before dN calculations
            obj.N0 = obj.N;

            % calculate dN from each enabled operator
            for ii = 1:length(operators)
                operatorToggle = char(operators(ii)); operatorToggle = [operatorToggle(3:end), 'Toggle'];
                if obj.(operatorToggle)
                    obj.N = obj.N + obj.(operators(ii));
                end
            end

            obj.N(obj.N >-1e-4 & obj.N<0) = 0; % remove slightly negative values

            %% check for stability
            if any(obj.T(:)<abs(dT(:))) || any(obj.mv(:) < abs(dm(:))) || any(obj.N<0,"all")
                obj.dt = obj.dt/2;
                error('Something went wrong, make sure dt is small enough and all state variables are set and in the correct units.')
            end

            %% update state variables
            if obj.waterVaporInteraction
                obj.T = obj.T + dT;
                obj.pv = (obj.mv - dm).*obj.Rv.*obj.T/100;
            end

            %% update the source number concentration (i.e., at cloud base.)
            if obj.nSteps < obj.st
                obj.N(obj.NP > 0) = obj.NP(obj.NP > 0);
            end

            %% update the number of steps taken
            obj.nSteps = obj.nSteps + 1;
        end

        function obj = exmiras(obj)
            %% constructor for exmiras class, inherits from thermo class
            % builds the droplet size distribution, terminal velocity, 
            % radar, and da classes

            % calculate dropsize distribution
            obj.De = logspace(log10(0.1),log10(obj.dMax),obj.nBins+1);
            % obj.De = linspace(0.1, obj.dMax, obj.nBins+1); % mm; edges of droplet size bins
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