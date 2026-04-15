function cooling = coolingConfig(propellant)
cooling.cp  = 4180;                                      % [J/kg/K]
cooling.rho = 997;                                       % [kg/m^3]
cooling.mu  = 0.89e-3;                                   % [Pa*s]
cooling.kWater   = 0.60;                                 % [W/m/K]
cooling.Pr  = cooling.cp * cooling.mu / cooling.kWater;  % Prandtl number

cooling.inletTemperature = 18 + 273.15;             % [K]
cooling.kWall = 10;                                 % [W/m*K]
cooling.pressure = 10*1e5;                          % [Pa]
cooling.velocity = 10;                              % [m/s]

% Defining a struct to be used as cooling.
% Having unique struct for conductive and convective heat transfer

%% Fluid
cooling(1).type = "convective";
cooling(1).cp  = propellant.cea.cp;
cooling(1).rho = propellant.cea.rhoGas;
cooling(1).mu  = propellant.cea.mu;
cooling(1).k   = propellant.cea.k;
cooling(1).Pr  = cooling(1).cp * cooling(1).mu / cooling(1).k;
cooling(1).dx  = nan;
cooling(1).temperature = nan;
cooling(1).pressure = nan;
cooling(1).mDot = nan;

%% TBC
cooling(2).type = "conductive";
cooling(2).cp  = 4180;
cooling(2).rho = 997;
cooling(2).mu  = 0.89e-3;
cooling(2).k   = 0.60;
cooling(2).Pr  = cooling(2).cp * cooling(2).mu / cooling(2).k;
cooling(2).dx  = nan;
cooling(2).temperature = nan;
cooling(2).pressure = nan;
cooling(2).mDot = nan;

%% Wall
cooling(3).type = "conductive";
cooling(3).cp  = 4180;
cooling(3).rho = 997;
cooling(3).mu  = 0.89e-3;
cooling(3).k   = 0.60;
cooling(3).Pr  = cooling(3).cp * cooling(3).mu / cooling(3).k;
cooling(3).dx  = nan;
cooling(3).temperature = nan;
cooling(3).pressure = nan;
cooling(3).mDot = nan;

%% Coolant (Water)
cooling(4).type = "convective";
cooling(4).cp  = 4180;
cooling(4).rho = 997;
cooling(4).mu  = 0.89e-3;
cooling(4).k   = 0.60;
cooling(4).Pr  = cooling(4).cp * cooling(4).mu / cooling(4).k;
cooling(4).dx  = nan;
cooling(4).temperature = 273.15 + 18;
cooling(4).pressure = 10*1e5;
cooling(4).mDot = 10;
end

