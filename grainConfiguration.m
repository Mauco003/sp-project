function [dExt, L0] = grainConfiguration(dInt, web)

dExt = dInt + 2*web;
L0 = 1/2 * (3*dExt + dInt);