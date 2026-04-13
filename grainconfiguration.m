function [dext, L0] = grainconfiguration(dint, web)

dext = dint + 2*web;
L0 = 1/2 * (3*dext + dint);