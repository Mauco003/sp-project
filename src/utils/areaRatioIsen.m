function ratio = areaRatioIsen(M, gamma)
% areaRatioIsen - Computes the isentropic area ratio A/A* for a given Mach number and specific heat ratio.
%
% SYNTAX:
%  ratio = areaRatioIsen(M, gamma)
%
% INPUT:
%  M - Mach number (scalar or vector)
%  gamma - Specific heat ratio (scalar)
%
% OUTPUT:
%  ratio - Isentropic area ratio A/A* corresponding to the input Mach number(s)
%
% NOTE:
%  Referred to CRITICAL CONDITIONS (denominator is throat area)

    ratio = (1 / M) * ...
        ((2/(gamma + 1))*(1 + (gamma - 1)/2* M^2))^((gamma + 1)/(2*(gamma - 1)));

end