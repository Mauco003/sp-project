function propellant = propellantConfig()
% Propellant configuration from empirical data and CEA

% Taulated data
% Pressure [Pa]
propellant.ccPressure = [10.1; 9.7; 10.3; 30.2; 31.0; 29.8; 50.2; 51.0; ...
            50.3; 70.2; 69.0; 69.8; 91.2; 89.1; 89.0]*1e5;

% Burning rate [m/s]
propellant.burnRate = [4.0; 3.8; 4.1; 5.6; 6.0; 5.7; 7.0; 7.2; ...
                7.1; 8.4; 8.3; 8.6; 8.8; 9.0; 9.2]*1e-3;

% from NASA CEA output 80% AP - 20% HTPB
propellant.cea.gamma = 1.2386*1.00032;                                            % Specific heat ratio [-]
propellant.cea.ccTemperature = 2342.92;                                           % Chamber temperature [K]
propellant.cea.molarMass = 22.087;
propellant.cea.mu = 0.78107*1e-4;       % [Pa*s]
propellant.cea.k  = 2.7189 *1e-1;       % [W/m*K] Steady state heat conductivity

propellant.cea.rhoAP = 1950;                                                               % Ammonium perchlorate density [kg/m^3]
propellant.cea.rhoHTPB = 913;                                                              % HTPB density [kg/m^3]
propellant.cea.rhoP = 1/(0.8/propellant.cea.rhoAP +  0.2/propellant.cea.rhoHTPB);                                        % Propellant density [kg/m^3]
end