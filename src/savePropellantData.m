% Empirical data
% Pressure [bar]
data.ccPressure = [10.1; 9.7; 10.3; 30.2; 31.0; 29.8; 50.2; 51.0; ...
            50.3; 70.2; 69.0; 69.8; 91.2; 89.1; 89.0];

% Burning rate [mm/s]
data.burnRate = [4.0; 3.8; 4.1; 5.6; 6.0; 5.7; 7.0; 7.2; ...
                7.1; 8.4; 8.3; 8.6; 8.8; 9.0; 9.2];

% from NASA CEA output 80% AP - 20% HTPB 
data.cea.gamma = 1.2386*1.00032;                                                     % Specific heat ratio [-]
data.cea.ccTemperature = 2342.92;                                                               % Chamber temperature [K]
data.cea.molarMass = 22.087;

save('./data/propellant.mat', '-struct', 'data');