function cf = cfIdeal(gamma, pe, pc)
% cfIdeal - Ideal thrust coefficient for OPTIMAL EXPANSION
%
% SYNTAX:
%  cf = cfIdeal(gamma, pe, pc)
% 
% INPUT:
%  gamma    - specific heat ratio [-]
%  pe       - exit pressure [Pa]
%  pc       - chamber pressure [Pa]
%
% OUTPUT:
%  cf       - thrust coefficient [-]

    term1 = (2*gamma.^2./(gamma - 1)) .* (2./(gamma + 1)).^((gamma + 1)./(gamma - 1));
    term2 = 1 - (pe./pc).^((gamma - 1)./gamma);
    cf = sqrt(term1.*term2);

end