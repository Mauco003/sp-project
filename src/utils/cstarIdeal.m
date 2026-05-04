function cstar = cstarIdeal(R, Tc, gamma)
% cstarIdeal - Ideal characteristic velocity
%
% SYNTAX:
%  cstar = cstarIdeal(R, Tc, gamma)
%
% INPUT:
%  R      - specific gas constant [J/(kg*K)]
%  Tc     - chamber temperature [K]
%  gamma  - specific heat ratio [-]
%
% OUTPUT:
%  cstar  - characteristic velocity [m/s]

    cstar = sqrt(R*Tc/gamma) * ((gamma + 1)/2)^((gamma + 1)/(2*(gamma - 1)));

end