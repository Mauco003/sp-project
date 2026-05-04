function hg = bartzCorrelation(pc, cstar, Dt, rCurvature, epsilon, mu, cp, Pr, Twall, Tc, Te, omega)
% bartzCorrelation - compute the convective heat transfer coefficient using the Bartz correlation
% SYNTAX:
%  hg = bartzCorrelation(pc, cstar, Dt, rCurvature, epsilon
% INPUT:
%  pc         - chamber total pressure [Pa]
%  cstar      - characteristic velocity [m/s]
%  Dt         - throat diameter [m]
%  rCurvature - throat radius of curvature [m]
%  epsilon    - local area ratio A/At
%  mu         - gas viscosity [Pa·s]
%  cp         - gas cp [J/kg/K]
%  Pr         - gas Prandtl number [-]
%  Twall      - wall temperature [K]
%  Tc         - chamber / stagnation temperature [K]
%  Te         - local static temperature at the evaluated section [K]
%  omega      - viscosity exponent in Bartz correction [-]
%
% OUTPUT:
%  hg - convective heat transfer coefficient [W/m^2/K]
%
% Note:
% (1 + (gamma-1)/2 * M^2) = Tc / Te

    theta = Tc ./ Te;   % = 1 + (gamma-1)/2 * M^2

    sigma = 1 ./ ( ...
        (0.5 .* (Twall ./ Tc) .* theta + 0.5) .^ (0.8 - omega/5) .* ...
        theta .^ (omega/5) );

    hg = 0.026 ./ Dt.^0.2 .* (Dt ./ rCurvature).^0.1 .* ...
        (mu.^0.2 .* cp ./ Pr.^0.6) .* ...
        (pc ./ cstar).^0.8 .* ...
        (1 ./ epsilon).^0.9 .* ...
        sigma;

end