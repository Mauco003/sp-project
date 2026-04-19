function [t, performance, performanceCEA, grain] = computePerformance(a, n, propellant, performanceNom, nozzle, grain, constants, options)

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

%% Simplified Performance
% Mcc = machFromAreaRatio(nozzle.Acc/nozzle.At, gamma, 'subsonic'); % Computing mach in combustion chamber
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

%% CEA Performance
N = length(t);

for i=1:N
    temp_exit = getExitCEA(propellant.wtAP,propellant.wtHTPB, pc(i)*1e-5, nozzle.epsilon);
    exit_data.pressure(i) = temp_exit.pressure;               %needs BAR
    exit_data.mach(i)     = temp_exit.mach;
    exit_data.cf(i)       = temp_exit.cf;
    exit_data.cstar(i)    = temp_exit.cstar;
end

cfMomId = exit_data.cf;

cfMom = cfMomId .* lambda .* options.etaF;
cfStatic = (exit_data.pressure'-constants.pAmb)./pc * nozzle.epsilon;
cf = cfMom' + cfStatic;

cstar = exit_data.cstar;

T = mDot .* cstar' .* cf;
Isp = (cstar .* cf) ./ constants.g0;

ve = cstar .* cfMom;

performanceCEA.t = t;
performanceCEA.thrust = T;
performanceCEA.Isp = Isp;
performanceCEA.mDot = mDot;
performanceCEA.ve = ve;
performanceCEA.pc = pc;
performanceCEA.pe = exit_data.pressure;
performanceCEA.Me = exit_data.mach;
performanceCEA.rb = rb;
performanceCEA.cstar = exit_data.cstar;
performanceCEA.ct = exit_data.cf;

grain.dInt = 2*rInt;
grain.Ab = Ab;
grain.L = L;
end