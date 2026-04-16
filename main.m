%% TODO
% - computePerformance: add isoentropic efficiency
% - montecarlo: integrate in the rest of the toolkit and bufix
% - cooling jacket: implement CEA data inside computation
%
% - dimensioning of the cooling jacket

%% Setup
clear; close all; clc

addpath(genpath("./src"))

plotFlag = false;

propellant = propellantConfig();
cooling = coolingConfig(propellant);

%% Configuraiton Data
% Generic constants
constants = Constants();

% Requirements / design point
input.thrust = 100000;                                                            % [N]
input.totaImpulse = 2.5e6;                                                       % [N s]
input.burningTime = input.totaImpulse/input.thrust;                                          % [s]
input.pcNominal = 70e5;                                                                  % [Pa]
input.peNominal = constants.pAmb;                                                        % [Pa] (nozzle optimal at sea level)

% Area ratios (wrt At) from which to which there is active cooling
coolARatios = [2, 2];                                            

%% Ballistic characterization
% Vieille's law
[a, aSigma, n, nSigma, R2] = uncertaintyVieille(propellant.ccPressure, propellant.burnRate);

rbNominal = a*input.pcNominal^n;                                                         % Burning rate with nominal chamber pressure [mm/s]


%% Ideal thermodynamics & mass sizing
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
[performanceNom, AtIdeal, AeIdeal] = performanceNomCalc(input.thrust, propellant.cea.ccTemperature, input.pcNominal, input.peNominal, propellant.cea.gamma, propellant.cea.molarMass, constants);
=======
[performanceNom, AtIdeal, AeIdeal] = idealPerformance(thrust, propellant.cea.ccTemperature, pcNominal, peNominal, propellant.cea.gamma, propellant.cea.molarMass, constants);
>>>>>>> c52e2b5 (Reorganize into folders and cleanup)
=======
[performanceNom, AtIdeal, AeIdeal] = idealPerformance(input.thrust, propellant.cea.ccTemperature, input.pcNominal, input.peNominal, propellant.cea.gamma, propellant.cea.molarMass, constants);
>>>>>>> 5011f48 (Working nozzleThermalModel in main)
=======
[performanceNom, AtIdeal, AeIdeal] = performanceNomCalc(input.thrust, propellant.cea.ccTemperature, input.pcNominal, input.peNominal, propellant.cea.gamma, propellant.cea.molarMass, constants);
>>>>>>> 26367ee (Adapt some functions to changes in dev)

mPTot = performanceNom.mDot*input.burningTime;                                          % Total propellant mass [kg]


%% BATES motor desing
% Grain sizing (using Richard Nakka formulation)
grain = grainConfiguration(rbNominal, input.burningTime, mPTot, propellant.cea);
if plotFlag, plotHollowCylinder(grain); end

% Nozzle geometry (conical) sizing
%
rcc = grain.dExt0 / 2;                                                       % Chamber radius from grain OD

% Conical nozzleDesign
nozzle = nozzleDesign(AtIdeal, AeIdeal, rcc, ...
                                    "alpha", 15 * pi/180, ...
                                    "beta", 30 * pi/180, "showSummary",true);


%% Combustion chamber / Internal Ballistics (Real Performance)
[t, performance, grain] = computePerformance(a, n, ...
    propellant, performanceNom, nozzle, grain, constants);

%% Nozzle cooling
best = nozzleThermalModel(input, propellant, nozzle, ...
    performance, cooling, constants,"makePlot",true);