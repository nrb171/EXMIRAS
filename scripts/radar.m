classdef radar < handle
    properties
        D % liquid droplet diameter in mm, corresponding to the size of the droplet at the bins
        De
        Dw
        M
        nBins = 250

        dpp % dual-pol preprocessing variables

        dMax = 10 % maximum diameter in mm
        lambda % radar wavelength in mm

        mu % shape parameter
        gamma % slope parameter
        N0 % intercept parameter


        epsilon = 62.1 - sqrt(-1)*32.0; % dielectric constant of liquid water
        Kw = 0.93; % reflectivity factor of water
    end

    %% radar methods
    methods
        function dpp = dualPolPreprocessing(radar, lambda)
            % various preprocessing steps to calculate dual-pol variables
            % keyboard
            r = 0.9951 + 0.02510*radar.D - 0.03644*radar.D.^2 + 0.005303*radar.D.^3 - 0.0002492*radar.D.^4;
            b = ((radar.D/1000).^3./r).^(1/3);
            a = b.*r;

            reflh = readmatrix('../data/LUTs/reflh.txt');
            reflv = readmatrix('../data/LUTs/reflv.txt');
            xsecth = readmatrix('../data/LUTs/xsecth.txt');
            xsectv = readmatrix('../data/LUTs/xsectv.txt');
            rhohva = readmatrix('../data/LUTs/rhohva.txt');
            rhohvb = readmatrix('../data/LUTs/rhohvb.txt');
            kdpl = readmatrix('../data/LUTs/kdp.txt');
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
        function Zhh = calcZhh(radar, N, dpp)
            % in: N in m^-3 mm^-1, D in mm
            % out: horizontal reflectivity in dBZ
            % keyboard
            N = reshape(N, [], radar.nBins);
        

            Zhh = 10*log10(...
                trapz(radar.D, dpp.reflh.*N, 2) ...
            );
           
            Zhh = Zhh;
           
        end

        function Zvv = calcZvv(radar, N, dpp)
            % in: N in m^-3 mm^-1, D in mm
            % out: vertical reflectivity in dBZ
            N = reshape(N, [], radar.nBins);
            % dpp = radar.dualPolPreprocessing(radar.lambda);
           
            Zvv = 10*log10(...
                trapz(radar.D, dpp.reflv.*N, 2) ...
            );
           
            Zvv = Zvv;
           
        end

        function rhohv = calcRhohv(radar, N, dpp)
            N = reshape(N, [], radar.nBins);
            
            a = dpp.rhohva;
            b = dpp.rhohvb;
            c = dpp.xsecth./(2*pi);
            d = dpp.xsectv./(2*pi);
            
            rho_hv = ... 
                sqrt(...
                    ((trapz(radar.D, N.*a,2)).^2 + ...
                    (trapz(radar.D, N.*b,2)).^2)) ./ ...
                    sqrt((trapz(radar.D, N.*c,2) .* trapz(radar.D, N.*d,2)) ...
            );
            rhohv = reshape(rho_hv, size(radar.N, [1:3]));
        end
            
        function kdp = calcKdp(radar, N, dpp)
            N = reshape(N, [], radar.nBins);
            kdp = trapz(radar.D, N.*dpp.kdpl,2);
            kdp = reshape(kdp, size(radar.N, [1:3]));
        end
       
    end
    
    %% helpers
    methods
        function N = reshapeN(radar,N)
            N = reshape(N, [], radar.nBins);
            
            N = mat2cell(N, ones([1,size(N, 1)]));
        end
    end
    %% other useful methods
    methods
         %% initialize N based on reflectivity and Zdr
        function [N,N0, mu, gamma] = estimateDSD(radar, Z, Zdr)
            % in: Z in dBZ, Zdr in dBz
            % out: N in m^-3 mm^-1
            % keyboard
            inds = find(~isinf(Z));
            mu = NaN(size(Z));
            gamma = NaN(size(Z));
            % keyboard
            for i = 1:length(inds)
                ind = inds(i);
                % keyboard
                fun = @(x) calcErr(radar, radar.dpp, Z(ind), Zdr(ind), initN(radar, radar.dpp, Z(ind), x(1), x(2)), x(1), x(2));
                
                [x,fv] =  fminsearchbnd(fun, [5, 2.5], [0, -1], [20, 5]);

                N(inds(i), :) = initN(radar, radar.dpp, Z(ind), x(1), x(2));
                mu(inds(i)) = x(2);
                gamma(inds(i)) = x(1);
            end
            radar.mu = mu;
            radar.gamma = gamma;
            N0 = getN0(radar, radar.dpp, Z, gamma, mu);
            radar.N0 = N0;

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

        function [lwp] = calcLWP(radar, N)
            % in: N in m^-3 mm^-1
            % out: kg/m^3
            N = radar.reshapeN(N);
            lwp = NaN(size(N));
            for i = 1:length(N)
                Ni = N{i};
                lwp(i) = trapz(radar.D, Ni.*radar.M, 2);
                
            end
        end
    end
    % constructor
    methods
        function radar = radar(lambda)
            % get dsd parameters
            radar.De = logspace(log10(0.1),log10(radar.dMax),radar.nBins+1);
            radar.D =  (radar.De(1:end-1) + radar.De(2:end))/2;
            radar.Dw = diff(radar.De);

            rhol = 995.65; % kg/m^3 density of liquid water
            radar.M = (4*pi/3*rhol*(radar.D/1000/2).^3);


            radar.lambda = convertLambda(lambda);

            radar.dpp = radar.dualPolPreprocessing(radar.lambda);
            

            function lambdaOut = convertLambda(lambdaIn)
                if isa(lambda, 'str') | isa(lambda, 'char')
                    bandName = ["S",    "C",    "X",        "Ku",       "Ka"       ];
                    wavelength = [111,  53.5,   33.3,       22,         8.43       ];
                    lambdaOut = wavelength(strcmp(bandName, lambdaIn));
                else
                    lambdaOut = lambdaIn;
                end
            end % end convertLambda
        end % end constructor
    end % end constructor methods
end % end Radar class