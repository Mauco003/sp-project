%% Setup
% Maybe is not a professional thing, but here a playlist we listened to while coding this is:
% https://open.spotify.com/playlist/7rv59VIUexwUDMbiUhJUmb?si=900d566736ff42dc
% Hope this playlist make the work less hard and more enjoyable for you as it did for us :)

clear; close all; clc

addpath(genpath("."))

propellant = propellantConfig();
cooling = coolingConfig(propellant);
[casing, liner] = casingConfig();

x

%% Configuraiton Data
% Generic constants
constants = Constants();
design = designConfig(constants);

input = design.input;
plotFlag = design.plot.enableFigures;
savePlotsFlag = design.plot.savePlots;

%% Ballistic characterization
% Vieille's law
[a, aSigma, n, nSigma, R2] = uncertaintyVieille(propellant.ccPressure, propellant.burnRate);

rbNominal = a*input.pcNominal^n;                                           % Burning rate with nominal chamber pressure [mm/s]


%% Ideal thermodynamics & mass sizing
[performanceNom, AtIdeal, AeIdeal] = computeNominalPerformance(input.thrust, propellant.cea.ccTemperature, input.pcNominal, input.peNominal, propellant.cea.gamma, propellant.cea.molarMass, constants);
mPTot = performanceNom.mDot * input.burningTime;                           % Total propellant mass [kg]

%% BATES motor desing
% Grain sizing (using Richard Nakka formulation)
grain = grainDesign(rbNominal, input.burningTime, mPTot, propellant.cea);

% Conical nozzleDesign
nozzle = nozzleDesign(AtIdeal, AeIdeal, grain.dExt0 / 2, ...
                                    "alpha", design.nozzle.alpha, ...
                                    "beta", design.nozzle.beta, "showSummary",true);

%% Combustion chamber / Internal Ballistics (Real Performance)
[t, performance, performanceCEA, grain] = computePerformance(a, n, ...
    propellant, performanceNom, nozzle, grain, constants);

%% Nozzle cooling

thermalDesign = thermalModelDesign(input, nozzle, performanceCEA, cooling, ...
    "savePlots", savePlotsFlag);
thermalModelOut = nozzleThermalModel(thermalDesign, input, propellant, nozzle, ...
    performanceCEA, cooling, constants, "showSummary", true, "savePlots", savePlotsFlag);

%% Combustion chamber sizing
cc = combustionChamberDesign(input, performanceNom, grain, nozzle, ...
    propellant, casing, liner, constants);

%% Plots
if plotFlag
    plotHollowCylinder(grain);
    plotPerformance(performance, performanceCEA, "savePlots", savePlotsFlag);
end