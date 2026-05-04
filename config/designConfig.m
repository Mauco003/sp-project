function design = designConfig(constants)
% designConfig - Configuration data for the design of a BATES motor
%
% SYNTAX:
%  design = designConfig(constants)
%
% INPUT:
%  constants - Constants object containing physical constants
%
% OUTPUT:
%  design - Struct containing design configuration data

    arguments
        constants Constants = Constants()
    end

    % Design configuration
    design.input.thrust = 100000;                                               % [N]
    design.input.totalImpulse = 2.5e6;                                          % [N s]
    design.input.burningTime = design.input.totalImpulse / design.input.thrust; % [s]
    design.input.pcNominal = 70e5;                                              % [Pa]
    design.input.peNominal = constants.pAmb;                                    % [Pa]

    % Nozzle design parameters
    design.nozzle.alpha = 15 * pi / 180;
    design.nozzle.beta = 30 * pi / 180;

    % Plotting configuration
    design.plot.enableFigures = true;
    design.plot.savePlots = false;
end