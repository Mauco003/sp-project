function [thrust, Isp, mDot] = computePerformance(a, n, data, performanceNom, nozzle, grain, constants)
    
    [t, pc, rb] = computeBurn(a, n, data.cea.rhoP, performanceNom.cstar, grain, nozzle.At);

    y = cumtrapz(t, rb);

    rInt = (grain.dInt/2)+y;
    L = grain.L0 - 2*y;


    Ab = 2*pi*(0.25*grain.dExt^2 - rInt.^2) + 2*pi*rInt.*L;

    mDot = Ab .* rb * data.cea.rhoP;
    
    % thrust = mDot .* performanceNom.cstar;
    R = constants.R / data.cea.molarMass;
    gamma = data.cea.gamma;
    tc = data.cea.ccTemperature;

    Mcc = machFromAreaRatio(nozzle.Acc/nozzle.At, gamma, 'subsonic'); % Computing mach in combustion chamber
    Me = machFromAreaRatio(nozzle.Ae/nozzle.At, gamma, 'supersonic'); % Computing mach at nozzle exit
    
    p0 = pc*(1+0.5*(gamma-1)*Mcc^2)^(gamma/(gamma-1)); % Total cc pressure

    pe = p0 / (1 + 0.5*(gamma-1)*Me^2)^(gamma/(gamma-1)); % Static pressure in combustion chamber

    ve = exhaustVelocityIdeal(gamma, R, tc, pe, pc);
    Isp = ve/constants.g0; % Specific impulse
    thrust = mDot .* ve + (pe - data.cea.atmosphericPressure) * nozzle.Ae; % Thrust
end