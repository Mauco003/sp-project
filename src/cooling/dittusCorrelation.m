function hw = dittusCorrelation(Dc,k,Re,Pr)

% Dc in meters
% k % [W/m*K]

Nu = 0.023*Re^0.8*Pr^0.4;

hw = Nu * k / Dc; 

end