% mainMC - Multiphysics Analysis of Grain regression, Gas properties and Internal cooling (MAGGI)
% A unified Monte Carlo simulation script combining comprehensive 
% performance analysis and rapid bounds analysis based on distinct modes.
% Includes an optional thermal risk analysis for the cooling jacket.

% CHANGE CONFIGURATION OF THE SIMULATION SOME SECTIONS BELOW


%% Setup

clear; close all; clc

% Run main.m first to build the environment
main;
tic;


%% Configuration of the simulation

% Select Uncertainty Mode:
% 0 = Global Simulation (Vary ALL parameters simultaneously)
% 1 = Baseline Ballistics (Vary a, n. Freeze geometry/chemistry)
% 2 = Manufacturing Geometric (Vary tolerances. Freeze ballistics/chemistry)
% 3 = Thermodynamic (Vary O/F. Freeze geometry/ballistics)
mcConfig.analysisMode = 1; 

% Toggle Physical Computation Depth:
% true  = Deep Dive (Thrust, Isp, Tornado plots, Trade-off matrix)
% false = Fast Bounds (MEOP, Burn Time, Convergence tracking)
mcConfig.useComprehensive = false; 

% Toggle Thermal Finale:
% true = Run nozzleThermalModel using distributions (works only on mode 1)
mcConfig.runCoolingJacket = true; 

% Sampling Settings
mcConfig.correlationType = 'correlated'; % 'correlated' or 'independent'
mcConfig.numIter = 100; % Base iterations (generates numIter^2 simulations)

% Plotting Aesthetic Config
mcConfig.plots.plotFigs   = true;
mcConfig.plots.plotSumm   = true;
mcConfig.plots.fontName   = 'Times New Roman';
mcConfig.plots.fontSize   = 18;
mcConfig.plots.techBlue   = [0, 0.4470, 0.7410];       
mcConfig.plots.techOrange = [0.8500, 0.3250, 0.0980];  
mcConfig.plots.techGreen  = [0.176, 0.416, 0.310];     
mcConfig.plots.techPurple = [0.4, 0.2, 0.6];
mcConfig.plots.lineWidth  = 1.5;


%% Initialization & global Uncertainty generation

% Nominal Base Values
aNominal      = a;
nNominal      = n;

% Setup a & n sampling
switch mcConfig.correlationType
    case 'independent'
        [~, aStd, ~, nStd, ~] = uncertaintyVieille(propellant.ccPressure * 1e-6, propellant.burnRate * 1e3);
        aStd = aStd * 1e-3/(10^(6*nNominal));

        aMCUnshuffled = aNominal + aStd * randn(mcConfig.numIter, 1); 
        nMCUnshuffled = nNominal + nStd * randn(mcConfig.numIter, 1); 

        [aMCVec, nMCVec] = meshgrid(aMCUnshuffled, nMCUnshuffled);
        aMCVec = aMCVec(:);
        nMCVec = nMCVec(:);
        totalRunsMC = numel(aMCVec);
    case 'correlated'
        [~, ~, ~, qmMu, qmCovMatrix] = uncertaintyVieilleDos(propellant.ccPressure, propellant.burnRate);
        totalRunsMC = mcConfig.numIter^2;
        qmSamples = mvnrnd(qmMu, qmCovMatrix, totalRunsMC);
        
        qMC = qmSamples(:, 1);
        nMCVec = qmSamples(:, 2);
        aMCVec = exp(qMC);
    otherwise
        error("Invalid correlationType.");
end


%% Global Uncertainty Generation

% Nominal Values
OFNominal     = 80/20;
pAmbNominal   = constants.pAmb;
dtNominal     = sqrt(4/pi * nozzle.At);
deNominal     = sqrt(4/pi * nozzle.Ae);
dExt0Nominal  = grain.dExt0;
dInt0Nominal  = grain.dInt0;
L0Nominal     = grain.L0;
betaNominal   = design.nozzle.beta; 
alphaNominal  = design.nozzle.alpha;

% Standard Deviations (3-sigma spread)
deltaWall    = 0.00001;
LWallDiv     = nozzle.x(3)/cos(alphaNominal);
OFStd    = (0.05 * OFNominal) / 3;     
pAmbStd  = (0.04 * pAmbNominal) / 3;   
dtStd    = (0.01 * dtNominal) / 3;     
deStd    = (0.01 * deNominal) / 3;     
dExtStd  = (0.005 * dExt0Nominal) / 3;  
dIntStd  = (0.005 * dInt0Nominal) / 3;  
LStd     = (0.005 * L0Nominal) / 3;
alphaStd = atan(deltaWall/LWallDiv) / 3;

% Pre-generate Random Vectors
OFVec    = OFNominal + randn(totalRunsMC, 1) * OFStd;
pAmbVec  = pAmbNominal + randn(totalRunsMC, 1) * pAmbStd;
AtVec    = (dtNominal + randn(totalRunsMC, 1) * dtStd).^2 * pi/4;
AeVec    = (deNominal + randn(totalRunsMC, 1) * deStd).^2 * pi/4;
dExtVec  = dExt0Nominal + randn(totalRunsMC, 1) * dExtStd;
dIntVec  = dInt0Nominal + randn(totalRunsMC, 1) * dIntStd;
LVec     = L0Nominal + randn(totalRunsMC, 1) * LStd;
alphaVec = alphaNominal + randn(totalRunsMC, 1) * alphaStd;
MOxVec   = OFVec ./ (1 + OFVec);
MFVec    = 1 ./ (1 + OFVec);


%% 2D CEA Surrogate model generation

% The CEA outputs valuable differences are caught only by OF and epsilon variation
% It is assumed that the variation in pressure are enough small to be considered
% negligible for a frozen evolution simulation by CEA
nPoints = 10;

OFRange  = linspace(min(OFVec)*0.9, max(OFVec)*1.1, nPoints);
OFRange  = sort([OFRange, OFNominal]);

epsVec   = AeVec ./ AtVec;
epsRange = linspace(min(epsVec)*0.9, max(epsVec)*1.1, nPoints);
epsRange = sort([epsRange, (deNominal/dtNominal)^2]);

nPoints = length(epsRange);
[OFGrid, epsGrid] = ndgrid(OFRange, epsRange);

surrGamma   = zeros(nPoints, nPoints); surrR     = zeros(nPoints, nPoints);
surrTcc     = zeros(nPoints, nPoints); surrCstar = zeros(nPoints, nPoints);
surrPeRatio = zeros(nPoints, nPoints); surrMach  = zeros(nPoints, nPoints);
surrCfMom   = zeros(nPoints, nPoints);

% Building the surrogates of CEA (smooth functions)
parfor i = 1:nPoints
    MOxtmp = OFRange(i) / (1 + OFRange(i)) * 100;
    MFtmp  = 100 - MOxtmp;
    for j = 1:nPoints
        [chamberDatatmp, gasDatatmp] = getThermoProfileCEAFroz(MOxtmp, MFtmp, input.pcNominal*1e-5, 2, epsRange(j));
        surrGamma(i,j)   = chamberDatatmp.gamma;
        surrR(i,j)       = 8314.46 / chamberDatatmp.molarMass; 
        surrTcc(i,j)     = chamberDatatmp.T;
        surrCstar(i,j)   = gasDatatmp.exit.cstar; 
        surrMach(i,j)    = gasDatatmp.exit.mach;
        surrPeRatio(i,j) = gasDatatmp.exit.pressure / (input.pcNominal * 1e5); 
        surrCfMom(i,j)   = gasDatatmp.exit.cf;
    end
end

getGamma   = griddedInterpolant(OFGrid, epsGrid, surrGamma, 'spline');
getR       = griddedInterpolant(OFGrid, epsGrid, surrR, 'spline');
getTcc     = griddedInterpolant(OFGrid, epsGrid, surrTcc, 'spline');
getCstar   = griddedInterpolant(OFGrid, epsGrid, surrCstar, 'spline');
getMach    = griddedInterpolant(OFGrid, epsGrid, surrMach, 'spline');
getPeRatio = griddedInterpolant(OFGrid, epsGrid, surrPeRatio, 'spline');
getCfMom   = griddedInterpolant(OFGrid, epsGrid, surrCfMom, 'spline');


%% Execution loop (parfor)

% Pre-allocate actual applied parameter trackers
aMCActual = zeros(totalRunsMC, 1); nMCActual = zeros(totalRunsMC, 1);
OFMCActual = zeros(totalRunsMC, 1); pAmbMCActual = zeros(totalRunsMC, 1);
AtMCActual = zeros(totalRunsMC, 1); AeMCActual = zeros(totalRunsMC, 1);
dExtMCActual = zeros(totalRunsMC, 1); dIntMCActual = zeros(totalRunsMC, 1);
LMCActual = zeros(totalRunsMC, 1); alphaMCActual = zeros(totalRunsMC, 1);

% Pre-allocate always-on outputs
MEOP = zeros(totalRunsMC, 1);
burnTime = zeros(totalRunsMC, 1);

% Pre-allocate other outputs
if mcConfig.useComprehensive
    TMax = zeros(totalRunsMC, 1); TAvg = zeros(totalRunsMC, 1);
    IspMax = zeros(totalRunsMC, 1); IspAvg = zeros(totalRunsMC, 1);
    IspTotMax = zeros(totalRunsMC, 1); thrustReliability = zeros(totalRunsMC, 1);
    rB = zeros(totalRunsMC, 1);
end


% Parallel for execution for speed up simulations time
parfor index = 1:totalRunsMC

    % Set Nominal Baseline (the parfor require variables independent between each iteration)
    aMC = aNominal; nMC = nNominal;
    AtMC = nozzle.At; AeMC = nozzle.Ae; dExtMC = grain.dExt0; 
    dIntMC = grain.dInt0; LMC = grain.L0; pAmbMC = constants.pAmb; alphaMC = design.nozzle.alpha;
    OFMC = OFNominal; MOxMC = OFNominal/(1+OFNominal); MFMC = 1/(1+OFNominal);

    % Apply perturbations based on mcConfig.analysisMode
    switch mcConfig.analysisMode
        case 0 % Global (All variables)
            aMC = aMCVec(index); nMC = nMCVec(index);
            AtMC = AtVec(index); AeMC = AeVec(index); dExtMC = dExtVec(index);
            dIntMC = dIntVec(index); LMC = LVec(index); pAmbMC = pAmbVec(index);
            alphaMC = alphaVec(index);
            OFMC = OFVec(index); MOxMC = MOxVec(index); MFMC = MFVec(index);
        case 1 % Ballistics only
            aMC = aMCVec(index); nMC = nMCVec(index);
        case 2 % Geometric tolerances (+ ambient pressure) only
            AtMC = AtVec(index); AeMC = AeVec(index); dExtMC = dExtVec(index);
            dIntMC = dIntVec(index); LMC = LVec(index); pAmbMC = pAmbVec(index);
            alphaMC = alphaVec(index);
        case 3 % Thermodynamic only
            OFMC = OFVec(index); MOxMC = MOxVec(index); MFMC = MFVec(index);
    end
    

    epsMC = AeMC / AtMC;
    
    % Interpolate properties
    gammaMC = getGamma(OFMC, epsMC); 
    RMC     = getR(OFMC, epsMC); 
    TccMC   = getTcc(OFMC, epsMC); 
    cStarMC = getCstar(OFMC, epsMC);
    rhoPMC  = 1/(MOxMC/propellant.cea.rhoAP + MFMC/propellant.cea.rhoHTPB);
    
    surrogateData = [];
    surrogateData.Me      = getMach(OFMC, epsMC);
    surrogateData.PeRatio = getPeRatio(OFMC, epsMC);
    surrogateData.cfMom   = getCfMom(OFMC, epsMC);
        

    % Track actual applied values
    aMCActual(index) = aMC; nMCActual(index) = nMC;
    OFMCActual(index) = OFMC; pAmbMCActual(index) = pAmbMC;
    AtMCActual(index) = AtMC; AeMCActual(index) = AeMC;
    dExtMCActual(index) = dExtMC; dIntMCActual(index) = dIntMC;
    LMCActual(index) = LMC; alphaMCActual(index) = alphaMC;
    
    % Local struct packing
    grainMC = grain; grainMC.dExt0 = dExtMC; grainMC.dInt0 = dIntMC; grainMC.L0 = LMC;
    nozzleMC = nozzle; nozzleMC.At = AtMC; nozzleMC.Ae = AeMC;
    constantsMC = constants; constantsMC.pAmb = pAmbMC;

    propellantMC = propellant;
    propellantMC.cea.gamma = gammaMC; propellantMC.cea.R = RMC; 
    propellantMC.cea.ccTemperature = TccMC; propellantMC.cea.rhoP = rhoPMC;
    
    perfNomMC = performanceNom; perfNomMC.cstar = cStarMC;


    % Execution
    if mcConfig.useComprehensive
        % Full performance computation
        [tMCOutput, ~, perfMCOutput , ~] = computePerformance(aMC, nMC, ...
            propellantMC, perfNomMC, nozzleMC, grainMC, constantsMC, ...
            "alpha", alphaMC, "CEASurrogate", surrogateData);
        
        % Extract from output valuable info
        MEOP(index) = max(perfMCOutput.pc);
        burnTime(index) = tMCOutput(end) - tMCOutput(1);
        TMax(index) = max(perfMCOutput.thrust);
        TAvg(index) = mean(perfMCOutput.thrust);
        IspAvg(index) = mean(perfMCOutput.Isp);
        IspMax(index) = max(perfMCOutput.Isp);
        IspTotMax(index) = trapz(tMCOutput, perfMCOutput.thrust);
        rB(index) = aMC * input.pcNominal^nMC;
        
        % Evaluation of for how much time the thrust is "enough" (e.g. 90% of the design one)
        thrustThresholdRatio = 0.9;
        threshold = thrustThresholdRatio * design.input.thrust; 
        
        % Create a high-resolution time vector (2000 points provides excellent precision)
        tInterp = linspace(tMCOutput(1), tMCOutput(end), 2000);
        
        % Interpolate the thrust curve using 'pchip' to preserve shape without overshooting
        thrustInterp = interp1(tMCOutput, perfMCOutput.thrust, tInterp, 'pchip');
        
        % Find indices where the interpolated thrust meets the threshold
        idx = thrustInterp >= threshold;
        
        if any(idx)
            timeVec = tInterp(idx);
            % Calculates the continuous span from the first crossing to the last crossing
            thrustReliability(index) = ((timeVec(end) - timeVec(1)) / burnTime(index)) * 100;
        else
            % Fallback in case the motor fails to ever reach 90% of design thrust
            thrustReliability(index) = 0; 
        end

    else
        % Reduced performance evaluation (just time and pressure trace)
        [tMC, pMC, ~] = computeBurn(aMC, nMC, rhoPMC, cStarMC, grainMC, AtMC);
        MEOP(index) = max(pMC);
        burnTime(index) = tMC(end) - tMC(1);
    end
end

% Basic stats
MEOPMean = mean(MEOP); MEOPStd = std(MEOP);
burnTimeMean = mean(burnTime); burnTimeStd = std(burnTime);


%% Deterministic One-At-A-Time (OAT) marginal sensitivity

if mcConfig.analysisMode == 0 || mcConfig.analysisMode == 1
    
    numOatPoints = 100;
    
    % Extract the actual statistical spread the MC used
    stdA = std(aMCVec);
    stdN = std(nMCVec);
    
    % Generate perfectly spaced deterministic sweeps spanning +/- 3 sigma
    sweepA = linspace(aNominal - 3*stdA, aNominal + 3*stdA, numOatPoints)';
    sweepN = linspace(nNominal - 3*stdN, nNominal + 3*stdN, numOatPoints)';
    
    % Calculate exact percentage deviations for the X-axis
    pctDevA = (sweepA - aNominal) / aNominal * 100;
    pctDevN = (sweepN - nNominal) / nNominal * 100;
    MEOPVaryA = zeros(numOatPoints, 1); BTVaryA = zeros(numOatPoints, 1);
    MEOPVaryN = zeros(numOatPoints, 1); BTVaryN = zeros(numOatPoints, 1);
    
    cStarNom = performanceNom.cstar; rhoPNom = propellant.cea.rhoP; AtNom = nozzle.At;
    
    % Execute deterministic OAT sweeps
    parfor i = 1:numOatPoints
        % Vary 'a', hold 'n' strictly nominal
        [tA, pA, ~] = computeBurn(sweepA(i), nNominal, rhoPNom, cStarNom, grain, AtNom);
        MEOPVaryA(i) = max(pA); 
        BTVaryA(i) = tA(end) - tA(1);
        
        % Vary 'n', hold 'a' strictly nominal
        [tN, pN, ~] = computeBurn(aNominal, sweepN(i), rhoPNom, cStarNom, grain, AtNom);
        MEOPVaryN(i) = max(pN); 
        BTVaryN(i) = tN(end) - tN(1);
    end
end


%% Cooling Jacket risk analysis

% A steady-state themal model for the cooling jacket is used
% In order to get a valuable risk analysis we use the computed MEOP mean
% and std. The themal model needs just one value of chamber pressure, so
% the MEOP is chosen for testing

boilingLimit = cooling(4).Tboil;
TBCMaxTemperature = cooling(2).Tmax; % maxmimum temperature for TBC
INCONELMaxTemperature = cooling(3).Tmax;  % maxmimum temperature for INCONEL

if mcConfig.runCoolingJacket && mcConfig.analysisMode == 1

    ccPressureMCCJVec = MEOP;
    

    % After several analysis a mass flow rate of 1 kg/s was chosen
    % So a specification in the option of nozzleThermalModel was added
    mDotWaterMC = 1;
    if length(mDotWaterMC) > 1
        error("Further plots are not going to work if mDotWater is not a scalar)")
    else
        ncaseMDot = 1; 
    end

    % Pre-allocate outputs
    waterOutletTemperatureListMC = NaN(totalRunsMC, ncaseMDot);
    maxTcoldListMC = NaN(totalRunsMC, ncaseMDot);
    maxQCEAListMC = NaN(totalRunsMC, ncaseMDot);
    TBCHotListMC = NaN(totalRunsMC, ncaseMDot);
    INCONELHotListMC = NaN(totalRunsMC, ncaseMDot);


    parfor index = 1:totalRunsMC
        inputMC = input; inputMC.pcNominal = ccPressureMCCJVec(index);
        
        outMC = nozzleThermalModel(thermalDesign, inputMC, propellant, nozzle, ...
            performanceCEA, cooling, constants, "showSummary", false, ...
            "savePlots", false, "mDotList", mDotWaterMC, "GasCEA", thermalModelOut.GasCEA);

        TCEA_all = outMC.TCEA_all{:};
        waterOutletTemperatureListMC(index) = outMC.waterOutletTemperatureList;
        maxTcoldListMC(index) = outMC.maxTcoldList;
        maxQCEAListMC(index) = outMC.maxQCEAList;
        TBCHotListMC(index) = max(TCEA_all(:, 2));
        INCONELHotListMC(index) = max(TCEA_all(:, 3));
    end
    maxTcold = maxTcoldListMC;
    PoF = (sum(maxTcold >= boilingLimit) / totalRunsMC) * 100;
end


%% Data packing & plotting delegation


% Pack nominal inputs
mcData.nominals.burningTime = input.burningTime;
mcData.nominals.pcNominal = input.pcNominal;
mcData.nominals.a = aNominal;
mcData.nominals.n = nNominal;
mcData.totalRunsMC = totalRunsMC;

% Pack true simulated nominals from the baseline execution in main.m
mcData.nominals.simMEOP = max(performance.pc);
mcData.nominals.simBurnTime = t(end) - t(1);


% Pack base outputs
mcData.MEOP = MEOP;
mcData.burnTime = burnTime;
mcData.MEOPMean = MEOPMean;
mcData.MEOPStd = MEOPStd;
mcData.burnTimeMean = burnTimeMean;
mcData.burnTimeStd = burnTimeStd;

% Pack actual applied vectors
mcData.inputs.a = aMCActual; 
mcData.inputs.n = nMCActual;
mcData.inputs.OF = OFMCActual; 
mcData.inputs.pAmb = pAmbMCActual;
mcData.inputs.At = AtMCActual; 
mcData.inputs.Ae = AeMCActual;
mcData.inputs.dExt = dExtMCActual; 
mcData.inputs.dInt = dIntMCActual;
mcData.inputs.L = LMCActual; 
mcData.inputs.alpha = alphaMCActual;

% Pack comprehensive outputs
if mcConfig.useComprehensive
    mcData.TMax = TMax;
    mcData.TAvg = TAvg;
    mcData.IspMax = IspMax;
    mcData.IspAvg = IspAvg;
    mcData.IspTotMax = IspTotMax;
    mcData.thrustReliability = thrustReliability;
    mcData.rB = rB;
end

% Pack OAT outputs
if mcConfig.analysisMode == 0 || mcConfig.analysisMode == 1
    mcData.oat.pctDevA = pctDevA;
    mcData.oat.pctDevN = pctDevN;
    mcData.oat.MEOPVaryA = MEOPVaryA;
    mcData.oat.MEOPVaryN = MEOPVaryN;
    mcData.oat.BTVaryA = BTVaryA;
    mcData.oat.BTVaryN = BTVaryN;
end

% Pack thermal outputs
if mcConfig.runCoolingJacket && mcConfig.analysisMode == 1
    mcData.thermal.mDotWater = mDotWaterMC;
    mcData.thermal.TOutlet = maxTcold;
    mcData.thermal.maxQ = maxQCEAListMC;
    mcData.thermal.TBCHot = TBCHotListMC;
    mcData.thermal.INCONELHot = INCONELHotListMC;
    mcData.thermal.PoF = PoF;
    mcData.thermal.boilingLimit = boilingLimit;
    mcData.thermal.TBCMaxLimit = TBCMaxTemperature; 
    mcData.thermal.INCONELMaxLimit = INCONELMaxTemperature;
end

% Call external plotting function
if mcConfig.plots.plotFigs, plotMCResults(mcData, mcConfig); end

% Call external console summary printout
if mcConfig.plots.plotSumm, printMCSummary(mcData, mcConfig); end

%%
tEndDD = toc;
fprintf('\nSimulation complete in %.2f seconds.\n', tEndDD);