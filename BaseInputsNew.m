function in = BaseInputsNew()
%BASEINPUTS Baseline case for the solid rocket motor design project
% All units are SI unless otherwise stated.
%
% This file only stores assumptions / baseline data.
% It should NOT perform sizing calculations.
%
% Usage:
%   in = BaseInputs();

    % Generic constants
    con = Constants();

    % Requirements / design point
    in.req.thrust = 100000;                                                 % [N]
    in.req.totalImpulse = 2.5e6;                                            % [N s]
    in.req.burningTime = in.req.totalImpulse/in.req.thrust;                 % [s]
    in.chamb.pc = 70e3;                                                     % [Pa]
    in.chamb.pe = con.pAmb;                                                 % [Pa] (nozzle optimal at sea level)
    
    % Area ratios (wrt At) from which to which there is active cooling
    in.req.coolARatios = [2, 2];                                            

    % Empirical data
    % Pressure [bar]
    in.data.pressureDataBar = [10.1; 9.7; 10.3; 30.2; 31.0; 29.8; 50.2; 51.0; ...
                50.3; 70.2; 69.0; 69.8; 91.2; 89.1; 89.0];
    
    % Burning rate [mm/s]
    in.data.burningRateDataMMS = [4.0; 3.8; 4.1; 5.6; 6.0; 5.7; 7.0; 7.2; ...
                    7.1; 8.4; 8.3; 8.6; 8.8; 9.0; 9.2];
    
    % from NASA CEA output 74% AP - 26% HTPB 
    in.prop.gamma = 1.2386*1.00032;                                         % Specific heat ratio [-]
    in.chamb.Tc = 2342.92;                                                  % Chamber temperature [K]
    

end