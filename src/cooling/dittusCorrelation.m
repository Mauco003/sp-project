function hw = dittusCorrelation(Dc,kH2O,Re,Pr)

Nu = 0.023*Re^0.8*Pr^0.4;

hg = Nu * kH2O / Dc; 
end