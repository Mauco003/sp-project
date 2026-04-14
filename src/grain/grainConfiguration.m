function grain = grainConfiguration(rbNominal, burningTime, mPTot, cea)

    web = rbNominal*burningTime;
    Vprop = mPTot/cea.rhoP;

    f = @(dInt) pi/8 * (3*(dInt + 2*web)^3 - 3*dInt^2*(dInt + 2*web) + dInt*(dInt + 2*web)^2 - dInt^3) - Vprop;
    dInt = fzero(f, web);

    dExt = dInt + 2*web;
    L0 = 1/2 * (3*dExt + dInt);



    grain.Vprop = Vprop;
    grain.web = web;
    grain.dInt0 = dInt;
    grain.dExt0 = dExt;
    grain.L0 = L0;
end