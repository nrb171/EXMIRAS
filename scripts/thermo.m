classdef thermo
    %% initial environmental thermodynamic properties
    properties
        T % temperature (of air volume) field in Kelvin
        p % pressure (of air volume) field in hPa
        rhoa    % density of air
    end

    %% constant thermodynamic properties
    properties (Access = protected)
        Rv = 461 %J/K/kg gas constant for water vapor in air
        Ra = 287.05 %gas constant for dry air
        rhol = 995.65 % kg/m^3 density of liquid water
        p0 = 1000 % hPa, reference pressure
        cp = 1.006e3 % J/kg K ; % J/kg/K specific heat capacity of air
        L %characteristic length scale of droplet
        U % characteristic velocity scale of droplet
    end

    properties (Dependent)
        FK      % heat conduction
        FD      % heat diffusion
        Lv      % latent enthalpy of vaporization
        K       % thermal conductivity of air
        Dv      % diffusivity of water vapor
        es      % saturation vapor pressure
        nu      % kinematic viscosity of air
        

        %% dimensionless numbers
        NSc     % Schmidt number
        NPr     % Prandtl number
        NRe     % Reynolds number

    end

    % derived thermodynamic properties
    methods 
        %% terminal velocity of droplets
        function L = get.L(obj)
            %% subclasses must define D, the diameter of the droplet(s)
            L = obj.D/1000;
            L = (4/3*pi*(L/2).^3)./(4*pi*(L/2).^2);
        end
        function U = get.U(obj)
            %% subclasses must define vt
            % note, this scheme assumes that the particles instantaneously 
            % reach their terminal velocity and do not accelerate over time 
            % in any dimension.
            U = obj.vt;
        end
        %% latent enthalpy of vaporization
        function Lv = get.Lv(obj)
            % in: T in Kelvin
            % out: Lv in J/kg
            % Source: Kumjian and Ryzhkov (2012) p. 1265
            gamma = 0.167+(3.67e-4)*obj.T;
            Lv = 2.499e6*(273.15./obj.T).^gamma;
        end
        %% thermal conductivity of air
        function K = get.K(obj)
            % in: T in Kelvin
            % out: K in J/m/s/K
            % Source: Kumjian and Ryzhkov (2012) p. 1265
            K = (0.441635 + 0.0071*obj.T)*1e-2;
        end
        %% diffusivity of water vapor
        function Dv = get.Dv(obj)
            % in: T in Kelvin
            % out: Dv in m^2/s
            % Source: Kumjian and Ryzhkov (2012) p. 1265
            Dv = 2.11e-5*(obj.T./273.15).^1.94.*(obj.p0./(obj.p));
        end
        %% saturation vapor pressure 
        function es = get.es(obj)
            % in: T in Kelvin
            % out: es in hPa
            % Source: Kumjian and Ryzhkov (2012) p. 1265
            es = 2.53e9*exp(-5420./(obj.T));
        end

        %% kinematic viscosity of air
        function nu = get.nu(obj)
            % in: T in Kelvin
            % out: nu in m^2/s
            % Source: Kumjian and Ryzhkov (2012) p. 1265
            etaa = (0.379565 + 0.0049*obj.T)*1e-4;
            nu = (etaa./obj.rhoa);
        end

        %% density of air
        function rhoa = get.rhoa(obj)
            % in: T in Kelvin
            % out: rho in kg/m^3
            % ideal gas law
            rhoa = 100*obj.p./(obj.Ra*obj.T);
        end

        %% dimensionless numbers
      
        %% Schmidt number
        function NSc = get.NSc(obj)
            % in: T in Kelvin
            % out: NSc, unitless
            NSc = obj.nu./obj.Dv;
        end

        %% Prandtl number
        function NPr = get.NPr(obj)
            % in: T in Kelvin
            % out: NPr, unitless
            
            NPr = obj.nu./(obj.K./(obj.rhoa.*obj.cp));
        end

        %% Reynolds number
        function NRe = get.NRe(obj)
            % in: T in Kelvin
            % out: NRe, unitless
            
            alpha = obj.L .* obj.U; 
            sz = size(obj.T);
            NRe = alpha ./ (obj.nu(:));
            NRe = reshape(NRe, [sz, numel(alpha)]);
            
        end

        %% heat conduction and diffusion for air
        %% heat conduction
        function FK = get.FK(obj)
            % in: T in Kelvin
            % out: FK in W/m/K
            % Source: Kumjian and Ryzhkov (2012) p. 1265
            fh = 0.78 + 0.308*obj.NPr.^(1/3) .* obj.NRe.^(1/2);
            FK = (obj.Lv./(obj.Rv.*obj.T) - 1).*obj.Lv.*obj.rhol./(fh.*obj.K.*obj.T);

        end

        %% heat diffusion
        function FD = get.FD(obj)
            % in: T in Kelvin
            % out: FD in 1
            % Source: Kumjian and Ryzhkov (2012) p. 1265
            fv = 0.78 + 0.308*obj.NSc.^(1/3) .* obj.NRe.^(1/2);
            FD = (obj.rhol*obj.Rv).*obj.T./(fv.*obj.Dv.*obj.es*100);
        end
    end
end