function cooling = coolingConfig(propellant)
% coolingConfig - Defines the cooling configuration for the nozzle design, including both
% convective and conductive heat transfer properties.
%
% SYNTAX:
%  cooling = coolingConfig(propellant)
%
% INPUT:
%  propellant - A struct containing the properties of the propellant, including CEA data
%
% OUTPUT:
%  cooling - A struct array containing the cooling configuration for both convective and conductive heat transfer



    % Defining a struct to be used as cooling.
    % Having unique struct for conductive and convective heat transfer

    % Safety factor for thermal design
    cooling.SF = 1.5;

    % Fluid properties (combustion gases)
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

    % TBC - selected material YCS
    cooling(2).type = "conductive";
    cooling(2).k   = 0.9;
    cooling(2).Tmax = 1800 + 273; % max of YCS [K]
    cooling(2).tmax = 600e-3; %[K]

    % Wall - selected material inconel

    cooling(3).type = "conductive";
    cooling(3).k   = 24.2;
    cooling(3).Tmax = 1000 + 273; % boiling of inconel [K]
    cooling(3).YS = 600e6; % max allowable stress

    % Coolant (Water)
    cooling(4).type = "convective";
    cooling(4).Tboil = 150 + 273;
    cooling(4).cp  = 4180;
    cooling(4).rho = 997;
    cooling(4).mu  = 0.89e-3;
    cooling(4).k   = 0.60;
    cooling(4).Pr  = cooling(4).cp * cooling(4).mu / cooling(4).k;
    cooling(4).dx  = 5*1e-4;
    cooling(4).Tini = 273.15 + 18;
    cooling(4).pressure = 10*1e5;

end

