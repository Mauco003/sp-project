%% WATER DUCT DSIGN
burnTime = 24;
coolingTime = burnTime + 4;     % [s]

% Area = total heat transfer area [m^2]
TH2Oinitial = 18 + 273.15;      % [K]
TH2OfinalGuess = 100 + 273.15;  % [K]
cH2O = 4180;                    % [J/(kg*K)]

Qdot = qdot * Area;             % [W]  if qdot = heat flux [W/m^2]
Q = Qdot * coolingTime;         % [J]

massH2O = Q / (cH2O * (TH2OfinalGuess - TH2Oinitial));   % [kg]
massFlowH2O = massH2O / coolingTime;                     % [kg/s]

%% Velocity


