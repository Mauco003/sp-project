function [performanceNom, At, Ae] = idealPerformance(thrust, tc, pc, pe, gamma, molarMass, constants, options)
%IDEALTHERMODYNAMICS Computes ideal thermodynamic performance and required throat/exit areas
%
arguments
    thrust
    tc
    pc
    pe
    gamma
    molarMass
    constants   Constants = Constants()
    options.showSummary = false
end

g0    = constants.g0;
R     = constants.R / molarMass;

% Ideal thermodynamic performance
ve    = exhaustVelocityIdeal(gamma, R, tc, pe, pc);
cstar = cstarIdeal(R, tc, gamma);
Isp   = ve / g0;
mdot  = thrust/(Isp*g0);

% Throat and exit areas
cf = cfIdeal(gamma, pe, pc);

At = thrust / (pc * cf);
epsilon = computeEpsilon(gamma, pe, pc);
Ae = epsilon * At;

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
    fprintf('ve (ideal)              : %.6f m/s\n', ve);
    fprintf('c* (ideal)              : %.6f m/s\n', cstar);
    fprintf('Isp (ideal)             : %.6f s\n', Isp);
    fprintf('mdot (ideal)            : %.6f kg/s\n', mdot);
    fprintf('cf (ideal)              : %.6f [-]\n', cf);
    fprintf('====================================================\n\n');
end

end