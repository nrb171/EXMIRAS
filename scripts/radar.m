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

        rngToggle = false; % toggle for random number generator

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
        function Zhh = calcZhh(radar, N)
            % in: N in m^-3 mm^-1 along the 
            % out: horizontal reflectivity in dBZ

            N = reshape(N, [], radar.nBins);
            Zhh = 10*log10(...
                trapz(radar.D, radar.dpp.reflh.*N, 2) ...
            );
        end

        function Zvv = calcZvv(radar, N)
            % in: N in m^-3 mm^-1, D in mm
            % out: vertical reflectivity in dBZ
            
            N = reshape(N, [], radar.nBins);
            Zvv = 10*log10(...
                trapz(radar.D, radar.dpp.reflv.*N, 2) ...
            );
        end

        function Zdr = calcZdr(radar, N)
            % in: N in m^-3 mm^-1, D in mm
            % out: differential reflectivity in dB

            N = reshape(N, [], radar.nBins);
            Zhh = radar.calcZhh(N);
            Zvv = radar.calcZvv(N);
            Zdr = Zhh - Zvv;
        end

        function rhohv = calcRhohv(radar, N)
            % in: N in m^-3 mm^-1, D in mm
            % out: copolar correlation coefficient (unitless)
            N = reshape(N, [], radar.nBins);
            
            a = radar.dpp.rhohva;
            b = radar.dpp.rhohvb;
            c = radar.dpp.xsecth./(2*pi);
            d = radar.dpp.xsectv./(2*pi);
            
            rhohv = ... 
                sqrt(...
                    ((trapz(radar.D, N.*a,2)).^2 + ...
                    (trapz(radar.D, N.*b,2)).^2)) ./ ...
                    sqrt((trapz(radar.D, N.*c,2) .* trapz(radar.D, N.*d,2)) ...
            );
            
        end
            
        function kdp = calcKdp(radar, N)
            % in: N in m^-3 mm^-1, D in mm
            % out: specific differential phase in deg/km
            N = reshape(N, [], radar.nBins);
            kdp = trapz(radar.D, N.*radar.dpp.kdpl,2);
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