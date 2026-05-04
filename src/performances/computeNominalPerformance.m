function [performanceNom, At, Ae] = computeNominalPerformance(thrust, tc, pc, pe, gamma, molarMass, constants, options)
% computeNominalPerformance - Computes nominal performance
%
% SYNTAX:
%  [performanceNom, At, Ae] = computeNominalPerformance(thrust, tc, pc, pe, gamma, molarMass, constants, options)
%
% INPUT:
%  thrust         - Desired thrust [N]
%  tc             - Combustion chamber temperature [K]
%  pc             - Combustion chamber pressure [Pa]
%  pe             - Nozzle exit pressure [Pa]
%  gamma          - Specific heat ratio [-]
%  molarMass      - Molar mass of the exhaust gases [kg/mol]
%  constants      - Constants object containing physical constants
%  options        - Structure containing simulation options:
%                   .showSummary - Flag to print performance summary (default: false)
%                   .alpha       - Nozzle contraction angle in radians (default: 15 degrees)
%                   .etaF        - Nozzle efficiency (default: 0.95)
%
% OUTPUT:
%  performanceNom - Structure containing nominal performance parameters:
%                   .ve    - Exhaust velocity [m/s]
%                   .cstar - Characteristic velocity [m/s]
%                   .Isp   - Specific impulse [s]
%                   .mDot  - Mass flow rate [kg/s]
%                   .cf    - Thrust coefficient [-]
%  At             - Throat area [m^2]
%  Ae             - Exit area [m^2]

    arguments
        thrust
        tc
        pc
        pe
        gamma
        molarMass
        constants   Constants = Constants()
        options.showSummary = false
        options.alpha = 15*pi/180
        options.etaF = 0.95
    end

    % Extract constants
    g0    = constants.g0;
    R     = constants.R / molarMass;

    % Compute the 2D correction factor for the conical nozzle design
    lambda = 0.5*(1+cos(options.alpha));


    % Ideal quantities
    cstarId = cstarIdeal(R, tc, gamma);
    cfId = cfIdeal(gamma, pe, pc);

    % Apply the efficiencies 
    cf = cfId*options.etaF*lambda;
    cstar = cstarId;

    % Compute the rest of the performance parameters
    Isp = cf*cstar/g0;
    mdot  = thrust/(cstar*cf);
    At = thrust / (pc * cf);
    epsilon = computeEpsilon(gamma, pe, pc);
    Ae = epsilon * At;
    ve = cstar*cf - (pe-constants.pAmb)/mdot*Ae;

    % Exporting data
    performanceNom.ve    = ve;
    performanceNom.cstar = cstar;
    performanceNom.Isp   = Isp;
    performanceNom.mDot  = mdot;
    performanceNom.cf    = cf;

    % Print summary
    if options.showSummary
        fprintf('\n');
        fprintf('====================================================\n');
        fprintf('              PERFORMANCE MODEL SUMMARY             \n');
        fprintf('====================================================\n');
        fprintf('ve (actual)             : %.6f m/s\n', ve);
        fprintf('c* (actual)             : %.6f m/s\n', cstar);
        fprintf('Isp (actual)            : %.6f s\n', Isp);
        fprintf('mdot (actual)           : %.6f kg/s\n', mdot);
        fprintf('cf (actual)             : %.6f [-]\n', cf);
        fprintf('====================================================\n\n');
    end

end