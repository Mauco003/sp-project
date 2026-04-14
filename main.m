%% Setup
clear; close all; clc

addpath(genpath("./src"))

plotFlag = false;

propellant = propellantConfig();
cooling = coolingConfig();

%% Configuraiton Data
% Generic constants
constants = Constants();

% Requirements / design point
thrust = 100000;                                                            % [N]
totalImpulse = 2.5e6;                                                       % [N s]
burningTime = totalImpulse/thrust;                                          % [s]
pcNominal = 70e5/1.035;                                                                  % [Pa]
peNominal = constants.pAmb;                                                        % [Pa] (nozzle optimal at sea level)

% Area ratios (wrt At) from which to which there is active cooling
coolARatios = [2, 2];                                            


%% Ballistic characterization
% Vieille's law
[a, aSigma, n, nSigma, R2] = uncertaintyVieille(propellant.ccPressure, propellant.burnRate);

rbNominal = a*pcNominal^n;                                                         % Burning rate with nominal chamber pressure [mm/s]


%% Ideal thermodynamics & mass sizing
<<<<<<< HEAD
[performanceNom, AtIdeal, AeIdeal] = performanceNomCalc(input.thrust, propellant.cea.ccTemperature, input.pcNominal, input.peNominal, propellant.cea.gamma, propellant.cea.molarMass, constants);
=======
[performanceNom, AtIdeal, AeIdeal] = idealPerformance(thrust, propellant.cea.ccTemperature, pcNominal, peNominal, propellant.cea.gamma, propellant.cea.molarMass, constants);
>>>>>>> c52e2b5 (Reorganize into folders and cleanup)

mPTot = performanceNom.mDot*burningTime;                                          % Total propellant mass [kg]


%% BATES motor desing
% Grain sizing (using Richard Nakka formulation)
grain = grainConfiguration(rbNominal, burningTime, mPTot, propellant.cea);
if plotFlag, plotHollowCylinder(grain); end



% Nozzle geometry (conical) sizing
%
rcc = grain.dExt0 / 2;                                                       % Chamber radius from grain OD

% Conical nozzleDesign
nozzle = nozzleDesign(AtIdeal, AeIdeal, rcc, ...
                                    "alpha", 15 * pi/180, ...
                                    "beta", 30 * pi/180);


%% Combustion chamber / Internal Ballistics (Real Performance)
[t, performance, grain] = computePerformance(a, n, ...
    propellant, performanceNom, nozzle, grain, constants);

%% Nozzle cooling

nozzleThermalModel(propellant, performance, nozzle, cooling, constants)