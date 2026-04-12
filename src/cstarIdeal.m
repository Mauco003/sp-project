function cstar = cstarIdeal(R, Tc, gamma)
% CSTARIDEAL
% Ideal characteristic velocity
%
% Inputs:
%   R      - specific gas constant [J/(kg*K)]
%   Tc     - chamber temperature [K]
%   gamma  - specific heat ratio [-]
%
% Output:
%   cstar  - characteristic velocity [m/s]

cstar = sqrt(R*Tc/gamma) * ((gamma + 1)/2)^((gamma + 1)/(2*(gamma - 1)));

end