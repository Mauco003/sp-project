function [casing, liner] = casingConfig()
% casingConfig - Configuration for the combustion chamber casing and liner
%
% sYNTAX:
%  [casing, liner] = casingConfig()
%
% INPUT
%  None
%
% OUTPUT
%  casing - Struct containing casing properties
%  liner - Struct containing liner properties

    % Casing properties
    casing.thermalConductivity = 42.7;  % Steel thermal conductivity
    casing.hoopStress = 460e6;          % update these w real vals
    casing.safetyFactor = 1.5;          % update w real values
    casing.TMax = 1432*0.5;
    casing.cost = 0.8;                  % [%/kg]
    casing.density = 7850;              % [kg/m^3]

    % Liner properties
    liner.thickness = 5e-3;
    liner.thermalConductivity = 0.225;
    liner.regressionRate = 0.00015;       % [kg/s*m^2]
    liner.cost = 7;                     % [$/kg]
    liner.density = 1208;               % [kg/m^3]
end

