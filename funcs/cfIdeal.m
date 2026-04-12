function cf = cfIdeal(gamma, pe, pc)
% CFIDEAL
% Ideal thrust coefficient for OPTIMAL EXPANSION
%
% Inputs:
%   gamma    - specific heat ratio [-]
%   pe       - exit pressure [Pa]
%   pc       - chamber pressure [Pa]
%
% Output:
%   CF       - thrust coefficient [-]

term1 = (2*gamma^2/(gamma - 1)) * (2/(gamma + 1))^((gamma + 1)/(gamma - 1));
term2 = 1 - (pe/pc)^((gamma - 1)/gamma);
cf = sqrt(term1*term2);

end