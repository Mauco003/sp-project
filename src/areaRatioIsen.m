function ratio = areaRatioIsen(M, gamma)
% AREARATIOISEN
% Isentropic area ratio A/A*
% Referred to CRITICAL CONDITIONS (denominator is throat area)

ratio = (1 / M) * ...
    ((2/(gamma + 1))*(1 + (gamma - 1)/2* M^2))^((gamma + 1)/(2*(gamma - 1)));

end