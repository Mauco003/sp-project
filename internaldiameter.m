function dint_num = internaldiameter(Vcc, web)

syms dint real
assume (dint>0)
dext = dint + 2*web;
eqn = pi/8 * (3*dext^3 - 3*dint^2*dext + dint*dext^2 - dint^3) == Vcc;
dint_num = solve(eqn, dint);
dint_num  = double(dint_num);