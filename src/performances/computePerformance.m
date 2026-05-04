function [t, performance, performanceCEA, grain] = computePerformance(a, n, propellant, performanceNom, nozzle, grain, constants, options)
% computePerformance -  Computes the performance history of a solid rocket motor based on the burn characteristics 
% and nozzle design, comparing a simplified model with CEA results.
%
% SYNTAX:
%  [t, performance, performanceCEA, grain] = computePerformance(a, n, propellant, performanceNom, nozzle, grain, constants, options)
%
% INPUT:
%  a                - Vieille's law coefficient [mm/s/bar^n]
%  n                - Vieille's law exponent
%  propellant       - A structure containing the propellant properties, including:
%                       .wtAP   - Weight percentage of Ammonium Perchlorate
%                       .wtHTPB - Weight percentage of Hydroxyl-terminated polybutadiene
%                       .cea    - A structure containing CEA properties of the propellant, including:
%                           .rhoP           - Propellant density [kg/m^3]
%                           .gamma          - Specific heat ratio [-]
%                           .R              - Specific gas constant [J/(kg*K)]
%                           .ccTemperature  - Combustion chamber temperature [K]
%  performanceNom   - A structure containing nominal performance parameters, including:
%                       .cstar - Nominal characteristic velocity [m/s]
%  nozzle           - A structure containing nozzle parameters, including:
%                       .At - Throat area [m^2]
%                       .Ae - Exit area [m^2]
%  grain            - A structure containing initial grain geometry, including:
%                       .dInt0 - Initial inner diameter [m]
%                       .dExt0 - Initial outer diameter [m]
%                       .L0 - Initial length [m]
%  constants        - A structure containing physical constants
%  options          - A structure containing simulation options:
%                       .alpha - Nozzle contraction angle in radians (default: 15 degrees)
%                       .etaF - Nozzle efficiency (default: 0.95)
%                       .etaTheta - Nozzle flow angular efficiency (default: 0.95)
%                       .CEASurrogate - A structure containing pre-computed surrogate data for CEA performance (default: [])
%
% OUTPUT:
%  t                - Time vector [s]
%  performance      - A structure containing the simplified performance parameters over time.
%  performanceCEA   - A structure containing the CEA performance parameters over time.
%  grain            - Updated grain geometry struct with fields:
%                       .dInt - Inner diameter history [m]
%                       .L    - Length history [m]

    arguments
        a
        n
        propellant
        performanceNom
        nozzle
        grain
        constants
        options.alpha = 15*pi/180
        options.etaF = 0.95
        options.etaTheta = 0.95
        options.CEASurrogate = []
    end

    % Compute the 2D correction factor for the conical nozzle design
    lambda = 0.5*(1+cos(options.alpha));

    % Compute the burn time and the pressure history
    [t, pc, rb, grain] = computeBurn(a, n, propellant.cea.rhoP, performanceNom.cstar, grain, nozzle.At);


    rInt = grain.rInt;
    L = grain.L;

    Ab = 2*pi*(0.25*grain.dExt0^2 - rInt.^2) + 2*pi*rInt.*L;
    Ab(Ab<0) = 0;

    mDot = Ab .* rb * propellant.cea.rhoP;


    % Simplified Performance

    gamma = propellant.cea.gamma;
    R = propellant.cea.R;
    Me = machFromAreaRatio(nozzle.Ae/nozzle.At, gamma, 'supersonic'); % Computing mach at nozzle exit

    p0 = pc; % Total cc pressure (Assumed to be same as static)
    pe = p0 ./ (1 + 0.5*(gamma-1)*Me.^2).^(gamma/(gamma-1)); % Static pressure in combustion chamber

    ve = exhaustVelocityIdeal(gamma, R, propellant.cea.ccTemperature, pe, pc);
    Isp = ve/constants.g0; % Specific impulse
    T = mDot .* ve + (pe - constants.pAmb) * nozzle.Ae; % Thrust

    performance.t = t;
    performance.thrust = T;
    performance.Isp = Isp;
    performance.mDot = mDot;
    performance.ve = ve;
    performance.pc = pc;
    performance.pe = pe;
    performance.Me = Me;
    performance.rb = rb;
    performance.cstar = performanceNom.cstar;
    performance.ct = performanceNom.cf;


    % CEA Performance
    pcbar = pc * 1e-5; % Convert pressure history to bar

    if isempty(options.CEASurrogate)
        % Standard mode: run localized CEA and create surrogate for next runs (if needed)
        nPoints = 10;
        pcMin = max(min(pcbar), 0.01); 
        pcMax = max(pcbar);
        pcSurrogate = linspace(pcMin, pcMax, nPoints);

        surrP     = zeros(nPoints, 1);
        surrMach  = zeros(nPoints, 1);
        surrCf    = zeros(nPoints, 1);

        for k = 1:nPoints
            tempExit = getExitCEA(propellant.wtAP, propellant.wtHTPB, pcSurrogate(k), nozzle.epsilon);
            surrP(k)     = tempExit.pressure;
            surrMach(k)  = tempExit.mach;
            surrCf(k)    = tempExit.cf;
        end

        exitData.pressure = interp1(pcSurrogate, surrP, pcbar, 'spline')';
        exitData.mach     = interp1(pcSurrogate, surrMach, pcbar, 'spline')';
        exitData.cf       = interp1(pcSurrogate, surrCf, pcbar, 'spline')';
    else
        % Fast mode: use 2D Pre-computed Surrogate Data
        peRatio  = options.CEASurrogate.PeRatio;
        exitMach = options.CEASurrogate.Me;
        cfMomId  = options.CEASurrogate.cfMom;
        
        % Expand constants dynamically over the pressure history
        exitData.pressure = pc .* peRatio;
        exitData.mach     = exitMach .* ones(size(pc));
        exitData.cf       = cfMomId .* ones(size(pc));
    end

    cfMomId = exitData.cf(:); % Force column vector for safety

    cfMom = cfMomId .* lambda .* options.etaF;
    cfStatic = (exitData.pressure(:) - constants.pAmb) ./ pc(:) * nozzle.epsilon;
    cf = cfMom + cfStatic;

    % Assumption of using nominal value of cstar
    cstar = performanceNom.cstar;

    T = mDot .* cstar .* cf;
    Isp = (cstar .* cf) ./ constants.g0;
    ve = cstar .* cfMom;


    % Pack CEA performance struct

    performanceCEA.t = t;
    performanceCEA.thrust = T;
    performanceCEA.Isp = Isp;
    performanceCEA.mDot = mDot;
    performanceCEA.ve = ve;
    performanceCEA.pc = pc;
    performanceCEA.pe = exitData.pressure;
    performanceCEA.Me = exitData.mach;
    performanceCEA.rb = rb;
    performanceCEA.cstar = performanceNom.cstar;
    performanceCEA.ct = exitData.cf;

    grain.dInt = 2*rInt;
    grain.Ab = Ab;
    grain.L = L;
end