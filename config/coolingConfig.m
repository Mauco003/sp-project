function cooling = coolingConfig(propellant)
cooling.kWall = 10;                                 % [W/m*K]

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
cooling(4).dx  = 5*1e-4;
cooling(4).temperature = 273.15 + 18;
cooling(4).pressure = 10*1e5;

cooling(4).dP = 0.5e5; % [Pa]
cooling(4).Dc   = 2*5e-3;                         % [m]
cooling(4).A = (cooling(4).Dc/2)^2*pi;               % [m^2]

cooling(4).velocity = sqrt(2*cooling(4).dP/cooling(4).rho);                           % [kg/s]
cooling(4).mDot = cooling(4).velocity*(cooling(4).A*cooling(4).rho); % [m/s]
cooling(4).Re = cooling(4).velocity*cooling(4).rho*cooling(4).Dc/cooling(4).mu; % [-]

end

