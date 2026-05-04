function hw = dittusCorrelation(Dc,k,Re,Pr)
% dittusCorrelation - Compute the Dittus correlation for convective heat transfer
%
% SYNTAX:
%  hw = dittusCorrelation(Dc,k,Re,Pr)
%
% INPUT:
%  Dc - Diameter of the tube [m]
%  k - Thermal conductivity of the fluid [W/m*K]
%  Re - Reynolds number
%  Pr - Prandtl number
% 
% OUTPUT:
%  hw - Convective heat transfer coefficient [W/m^2*K]

    Nu = 0.023*Re^0.8*Pr^0.4;

    hw = Nu * k / Dc; 

end