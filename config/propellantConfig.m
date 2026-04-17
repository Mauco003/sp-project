function propellant = propellantConfig(constants)
    arguments
        constants Constants = Constants()
    end
% Propellant configuration from empirical data and CEA

% Taulated data
% Pressure [Pa]
propellant.ccPressure = [10.1; 9.7; 10.3; 30.2; 31.0; 29.8; 50.2; 51.0; ...
            50.3; 70.2; 69.0; 69.8; 91.2; 89.1; 89.0]*1e5;

% Burning rate [m/s]
propellant.burnRate = [4.0; 3.8; 4.1; 5.6; 6.0; 5.7; 7.0; 7.2; ...
                7.1; 8.4; 8.3; 8.6; 8.8; 9.0; 9.2]*1e-3;

% from NASA CEA output 80% AP - 20% HTPB
wtAP = 80;
wtHPTB = 20;
p_chamber = 70; % bar
  eps_inlet = [2.0];
  eps_exit  = [2.0];

  
  [Chamber, ~] = getThermoProfileCEA_froz(wtAP, wtHPTB, p_chamber, eps_inlet, eps_exit);
  
propellant.cea.gamma = Chamber.gamma;                                            % Specific heat ratio [-]
propellant.cea.ccTemperature = Chamber.T;                                           % Chamber temperature [K]
propellant.cea.molarMass = Chamber.molar_mass;
propellant.cea.cp = Chamber.cp;                                                   % [J/kg/K]
propellant.cea.mu = Chamber.viscosity;       % [Pa*s]
propellant.cea.k  = Chamber.conductivity;       % [W/m*K] Steady state heat conductivity
propellant.cea.rhoGas = Chamber.density;     % [kg/m^3] Gas density at chamber conditions

propellant.cea.rhoAP = 1950;         % DATA -NOT- FROM CEA                                                     % Ammonium perchlorate density [kg/m^3]
propellant.cea.rhoHTPB = 913;                                                              % HTPB density [kg/m^3]
propellant.cea.rhoP = 100/(wtAP/propellant.cea.rhoAP +  wtHPTB/propellant.cea.rhoHTPB);                                        % Propellant density [kg/m^3]
end