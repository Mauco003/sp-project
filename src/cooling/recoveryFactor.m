function r = recoveryFactor(gamma, mach, Pr)
r = (1 + 0.5*Pr^(1/3)*(gamma-1)*mach.^2)./ ...
    (1+ 0.5*(gamma-1)*mach.^2);

end