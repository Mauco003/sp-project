function grain = grainDesign(rbNominal, burningTime, mPTot, cea)
% grainDesign - Computes the grain geometry for a BATES motor design
%
% SYNTAX:
%  grain = grainDesign(rbNominal, burningTime, mPTot, cea)
%
% INPUT:
%  rbNominal    - Nominal burning rate [mm/s]
%  burningTime  - Total burning time [s]
%  mPTot        - Total propellant mass [kg]
%  cea          - CEA data structure containing propellant properties
%
% OUTPUT:
%  grain        - Structure containing grain geometry parameters
    
    % Hp: the web is consumed at a nominal rate
    web = rbNominal*burningTime;
    Vprop = mPTot/cea.rhoP;

    % Solving for the internal diameter that satisfies the propellant volume requirement
    f = @(dInt) pi/8 * (3*(dInt + 2*web)^3 - 3*dInt^2*(dInt + 2*web) + dInt*(dInt + 2*web)^2 - dInt^3) - Vprop;
    dInt = fzero(f, web);

    dExt = dInt + 2*web;
    L0 = 1/2 * (3*dExt + dInt);


    % Pack output
    grain.Vprop = Vprop;
    grain.web = web;
    grain.dInt0 = dInt;
    grain.dExt0 = dExt;
    grain.L0 = L0;
end