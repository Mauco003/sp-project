function dintNum = internalDiameter(Vprop, web)

syms dInt real
assume (dInt>0)
dExt = dInt + 2*web;
eqn = pi/8 * (3*dExt^3 - 3*dInt^2*dExt + dInt*dExt^2 - dInt^3) == Vprop;
dintNum = solve(eqn, dInt);
dintNum  = double(dintNum);