function ve = exhaustVelocityIdeal(gamma, R, Tc, pe, pc)
% exhaustVelocityIdeal - Ideal exhaust velocity from isentropic expansion
%
% SYNTAX:
%  ve = exhaustVelocityIdeal(gamma, R, Tc, pe, pc)
%
% INPUT:
%  gamma  - specific heat ratio [-]
%  R      - specific gas constant [J/(kg*K)]
%  Tc     - chamber temperature [K]
%  pe     - exit pressure [Pa]
%  pc     - chamber pressure [Pa]
%
% OUTPUT:
%  ve     - ideal exhaust velocity [m/s]

    ve = sqrt((2*gamma./(gamma - 1)) .* R .* Tc .* (1 - (pe./pc).^((gamma - 1)./gamma)));

end