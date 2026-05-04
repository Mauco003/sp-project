classdef Constants
% Constants - Class to store constants used in the code
%
% SYNTAX:
%  constants = Constants()
%
% DESCRIPTION:
%  This class stores constants used in the code, such as the gravity constant, universal gas constant, and ambient temperature.
    
    properties (Constant)
        g0 = 9.81                                                           % Gravity constant [m/s^2]
        R  = 8314.29                                                        % Universal gas constant [J/kg/K]
        rb = 8.88;                                                           % Burning rate at nominal conditions [mm/s]
        tBurn = 25;                                                          % Burn time [s]
        TMax = 1000;                                                         % Maximum temperature [K]
        h_out = 10;                              % [W/m^2/K], natural convection
        T_ambient = 293.15;                      % [K], assumed
    end

    properties
        pAmb = 101325                                                       % Ambient pressure at sea level [Pa]
    end
end

