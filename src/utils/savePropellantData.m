function data = savePropellantData()
    % Script to save input data to .mat file

    % Empirical data
    % Pressure [Pa]
    data.ccPressure = [10.1; 9.7; 10.3; 30.2; 31.0; 29.8; 50.2; 51.0; ...
                50.3; 70.2; 69.0; 69.8; 91.2; 89.1; 89.0]*1e5;

    % Burning rate [m/s]
    data.burnRate = [4.0; 3.8; 4.1; 5.6; 6.0; 5.7; 7.0; 7.2; ...
                    7.1; 8.4; 8.3; 8.6; 8.8; 9.0; 9.2]*1e-3;

    % from NASA CEA output 80% AP - 20% HTPB
    data.cea.gamma = 1.2386*1.00032;                                            % Specific heat ratio [-]
    data.cea.ccTemperature = 2342.92;                                           % Chamber temperature [K]
    data.cea.molarMass = 22.087;
    data.cea.rhoAP = 1950;                                                               % Ammonium perchlorate density [kg/m^3]
    data.cea.rhoHTPB = 913;                                                              % HTPB density [kg/m^3]
    data.cea.rhoP = 1/(0.8/data.cea.rhoAP +  0.2/data.cea.rhoHTPB);                                        % Propellant density [kg/m^3]

    save('./data/propellant.mat', '-struct', 'data');

end