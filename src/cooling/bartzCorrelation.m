function hg = bartzCorrelation(pc, cstar, Dt, rCurvature, epsilon, mu, cp, Pr)

% Pc - TOtal pressure (chamber)
% rct - Radius of curvature of throat
% Dt - throat diameter
% Pr - Prandlt CEA
% epsilon - Area ratio at that point
% cstar
% cp from CEA
% mu from CEA (stagnation)

% densRef = (densE + densWall)/2;
% muRef = (muE + muWall)/2;

% sigma = (densRef/densE)^(0.8)*(muRef/mu0)^(0.8);
sigma = 1;

hg = 0.026 / Dt^0.2 * (Dt/rCurvature)^0.1* ...
     (mu^0.2 * cp / Pr^0.6) * ...
     (pc ./ cstar).^0.8 .* ...
     (1/epsilon).^0.9 .* ...
     sigma;
end