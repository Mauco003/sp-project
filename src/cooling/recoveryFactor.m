function r = recoveryFactor(gamma, mach, Pr)
% recoveryFactor - compute the recovery factor for a given flow condition
%
% SYNTAX:
%  r = recoveryFactor(gamma, mach, Pr)
%
% INPUT:
%  gamma - specific heat ratio of the gas [-]
%  mach  - local Mach number [-]
%  Pr    - local Prandtl number [-]
%
% OUTPUT:
%  r - recovery factor [-]

    r = (1 + 0.5*Pr^(1/3)*(gamma-1)*mach.^2)./ ...
        (1+ 0.5*(gamma-1)*mach.^2);

end