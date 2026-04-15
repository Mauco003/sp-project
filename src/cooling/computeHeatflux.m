function [q, Tout] = computeHeatFlux(Tin, gamma, mach, PrGas, PrWater, dx, dA, ...
    input, nozzle, performance, cooling)
% Compute heat flux and temperature distribution in the wall and coolant
% using the Bartz correlation for the gas-side heat transfer coefficient and

Tout(1) = T0*recoveryFactor(gamma, mach, PrGas); % Adiabatic wall temperature at station i

h1 = bartzCorrelation(input.pcNominal, performance.cstar, 2*rt, nozzle.rCurvature, epsilon(i), muGas, cpGas, PrGas);
k2 = cooling(2).k;
k3 = cooling(3).k;
h4 = 4000; %dittusCorrelation(Dc,kH2O,Re,Pr);

H = 1/(1/h1 + dx(1)/k2 + dx(2)/k3 + 1/h4);

q = H * (Tin(1) - Tin(end));

% Update station tempeartures
Tin(2) = Tin(1) - q/h1;
Tin(3) = Tin(2) - q*dx(1)/k2;
Tin(4) = Tin(3) - q*dx(2)/k3;

% Update water temperature at station i+1
Tout(end) = Tin(end) + q*dA/(cooling(end).cp*cooling(end).mDot);
end