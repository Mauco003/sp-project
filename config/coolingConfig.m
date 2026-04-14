function cooling = coolingConfig()
cooling.cp  = 4180;                                      % [J/kg/K]
cooling.rho = 997;                                       % [kg/m^3]
cooling.mu  = 0.89e-3;                                   % [Pa*s]
cooling.kWater   = 0.60;                                 % [W/m/K]
cooling.Pr  = cooling.cp * cooling.mu / cooling.kWater;  % Prandtl number

cooling.inletTemperature = 18 + 273.15;             % [K]
cooling.kWall = 10;                                 % [W/m*K]
cooling.pressure = 10*1e5;                          % [Pa]
cooling.velocity = 10;                              % [m/s]

end

