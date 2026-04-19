function [t, performance, grain] = computePerformance(a, n, propellant, performanceNom, nozzle, grain, constants, options)
    
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
    end

    lambda = 0.5*(1+cos(options.alpha));


    [t, pc, rb] = computeBurn(a, n, propellant.cea.rhoP, performanceNom.cstar, grain, nozzle.At);

    y = cumtrapz(t, rb);

    rInt = (grain.dInt0/2)+y;
    L = grain.L0 - 2*y;


    Ab = 2*pi*(0.25*grain.dExt0^2 - rInt.^2) + 2*pi*rInt.*L;
    Ab(Ab<0) = 0;

    mDot = Ab .* rb * propellant.cea.rhoP;
    
    N = length(t);

    for i=1:N
    temp_exit = getExitCEA(propellant.wtAP,propellant.wtHTPB, pc(i), nozzle.epsilon);
    exit_data.pressure(i) = temp_exit.pressure;
    exit_data.mach(i)     = temp_exit.mach;
    exit_data.cf(i)       = temp_exit.cf;
    exit_data.cstar(i)    = temp_exit.cstar;
    end
  

    cfMomId = exit_data.cf;

    cfMom = cfMomId .* lambda .* options.etaF;
    cfStatic = (exit_data.pressure'-constants.pAmb*eye(N,1))./pc * nozzle.epsilon;
    cf = cfMom + cfStatic;

    cstar = exit_data.cstar';

    thrust = mDot .* cstar .* cf;
    Isp = (cstar .* cf) ./ constants.g0;

    ve = cstar .* cfMom;
    
    performance.t = t;
    performance.thrust = thrust;
    performance.Isp = Isp;
    performance.mDot = mDot;
    performance.ve = ve;
    performance.pc = pc;
    performance.pe = exit_data.pressure;
    performance.Me = exit_data.mach;
    performance.rb = rb;
    performance.cstar = exit_data.cstar;
    performance.ct = exit_data.cf;
    
    grain.dInt = 2*rInt;
    grain.Ab = Ab;
    grain.L = L;
end