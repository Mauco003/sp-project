function [t, performance, grain] = computePerformance(a, n, data, performanceNom, nozzle, grain, constants, options)
    
    arguments
        a
        n
        data
        performanceNom
        nozzle
        grain
        constants
        options.alpha = 15*pi/180
        options.etaF = 0.95
        options.etaTheta = 0.95
    end

    lambda = 0.5*(1+cos(options.alpha));


    [t, pc, rb] = computeBurn(a, n, data.cea.rhoP, performanceNom.cstar, grain, nozzle.At);

    y = cumtrapz(t, rb);

    rInt = (grain.dInt0/2)+y;
    L = grain.L0 - 2*y;


    Ab = 2*pi*(0.25*grain.dExt0^2 - rInt.^2) + 2*pi*rInt.*L;
    Ab(Ab<0) = 0;

    mDot = Ab .* rb * data.cea.rhoP;

    gamma = data.cea.gamma;

    Mcc = machFromAreaRatio(nozzle.Acc/nozzle.At, gamma, 'subsonic'); % Computing mach in combustion chamber
    Me = machFromAreaRatio(nozzle.Ae/nozzle.At, gamma, 'supersonic'); % Computing mach at nozzle exit
    
    p0 = pc .* (1+0.5*(gamma-1)*Mcc.^2).^(gamma/(gamma-1)); % Total cc pressure

    pe = p0 ./ (1 + 0.5*(gamma-1)*Me.^2).^(gamma/(gamma-1)); % Static pressure in combustion chamber

    cfMomId = cfIdeal(gamma, pe, pc);
    cfMom = cfMomId .* lambda .* options.etaF;
    cfStatic = (pe-constants.pAmb)./pc * nozzle.epsilon;
    cf = cfMom + cfStatic;

    cstar = performanceNom.cstar;
    thrust = mDot .* cstar .* cf;
    Isp = (cstar .* cf) ./ constants.g0;

    ve = cstar .* cfMom;
    
    performance.t = t;
    performance.thrust = thrust;
    performance.Isp = Isp;
    performance.mDot = mDot;
    performance.ve = ve;
    performance.pc = pc;
    performance.pe = pe;
    performance.Mcc = Mcc;
    performance.Me = Me;
    performance.rb = rb;
    performance.cstar = performanceNom.cstar;
    
    grain.dInt = 2*rInt;
    grain.Ab = Ab;
    grain.L = L;
end