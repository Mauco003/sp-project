function eps = computeEpsilon(gamma, pe, pc)
% EPSILON
% Expansion ratio
%
% Inputs:
%   gamma    - specific heat ratio [-]
%   pe       - exit pressure [Pa]
%   pc       - chamber pressure [Pa]
%
% Output:
%   epsilon       - expansion ratio [-]

term1 = ((gamma + 1)/2)^(1/(gamma - 1));
term2 = (pe/pc)^(1/gamma);
term3 = (gamma + 1)/(gamma - 1)*(1 - (pe/pc)^((gamma - 1)/gamma));

eps = 1/(term1*term2*sqrt(term3));

end