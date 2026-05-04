function propellant = propellantConfig(constants)
% propellantConfig - set configuration data about the propellant (data from CEA)
%
% SYNTAX:
%  propellant = propellantConfig(constants)
%
% INPUT:
%  - constants [1x1 Constants] (optional)
%
% OUTPUT:
%  propellant  [1x1]  Struct with fields:
%                       - ccPressure [15x1]
%                       - burnRate   [15x1]
%                       - wtAP       [15x1]
%                       - wtHTPB     [15x1]
%                       - pChamber   [15x1]
%                       - cea        [1x1]  Struct with fields:
%                                           - gamma
%                                           - ccTemperature
%                                           - molarMass
%                                           - R
%                                           - cp
%                                           - mu
%                                           - k
%                                           - rhoGas
%                                           - cstar
%                                           - rhoAP
%                                           - rhoHTPB
%                                           - rhoP
    arguments
        constants Constants = Constants()
    end


    % Tabulated data

    % Pressure [Pa]
    propellant.ccPressure = [10.1; 9.7; 10.3; 30.2; 31.0; 29.8; 50.2; 51.0; ...
                50.3; 70.2; 69.0; 69.8; 91.2; 89.1; 89.0]*1e5;

    % Burning rate [m/s]
    propellant.burnRate = [4.0; 3.8; 4.1; 5.6; 6.0; 5.7; 7.0; 7.2; ...
                    7.1; 8.4; 8.3; 8.6; 8.8; 9.0; 9.2]*1e-3;


    % From NASA CEA output 80% AP - 20% HTPB
    propellant.wtAP = 80;
    propellant.wtHTPB = 20;
    propellant.pChamber = 70; % bar
    epsInlet = [2.0];
    epsExit  = [2.0];

    % Extract CEA data
    [chamber, gasInfo] = getThermoProfileCEAFroz(propellant.wtAP, propellant.wtHTPB, propellant.pChamber, epsInlet, epsExit);
    

    % Pack
    propellant.cea.gamma = chamber.gamma;                                           % Specific heat ratio [-]
    propellant.cea.ccTemperature = chamber.T;                                       % chamber temperature [K]
    propellant.cea.molarMass = chamber.molarMass;
    propellant.cea.R = constants.R / chamber.molarMass;
    propellant.cea.cp = chamber.cp;                                                 % [J/kg/K]
    propellant.cea.mu = chamber.viscosity;                                          % [Pa*s]
    propellant.cea.k  = chamber.conductivity;                                       % [W/m*K] Steady state heat conductivity
    propellant.cea.rhoGas = chamber.density;                                        % [kg/m^3] Gas density at chamber conditions

    propellant.cea.cstar = gasInfo.exit.cstar;


    propellant.cea.rhoAP = 1950;                                                    % DATA -NOT- FROM CEA                                                     % Ammonium perchlorate density [kg/m^3]
    propellant.cea.rhoHTPB = 913;                                                   % HTPB density [kg/m^3]
    propellant.cea.rhoP = 100/(propellant.wtAP/propellant.cea.rhoAP +  ...
                                propellant.wtHTPB/propellant.cea.rhoHTPB);          % Propellant density [kg/m^3]
end