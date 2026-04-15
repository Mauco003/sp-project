function ve = exhaustVelocityIdeal(gamma, R, Tc, pe, pc)
% EXHAUSTVELOCITYIDEAL
% Ideal exhaust velocity from isentropic expansion
%
% Inputs:
%   gamma  - specific heat ratio [-]
%   R      - specific gas constant [J/(kg*K)]
%   Tc     - chamber temperature [K]
%   pe     - exit pressure [Pa]
%   pc     - chamber pressure [Pa]
%
% Output:
%   ve     - ideal exhaust velocity [m/s]

ve = sqrt((2*gamma./(gamma - 1)) .* R .* Tc .* (1 - (pe./pc).^((gamma - 1)./gamma)));

end