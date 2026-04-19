%% TODO
% - computePerformance: add isoentropic efficiency
% - montecarlo: integrate in the rest of the toolkit and bufix
% - cooling jacket: implement CEA data inside computation
%
% - dimensioning of the cooling jacket

%% Setup
clear; close all; clc

addpath(genpath("."))

plotFlag = false;

propellant = propellantConfig();
cooling = coolingConfig(propellant);
[casing, liner] = casingConfig();

%% Configuraiton Data
% Generic constants
constants = Constants();

% Requirements / design point
input.thrust = 100000;                                                     % [N]
input.totaImpulse = 2.5e6;                                                 % [N s]
input.burningTime = input.totaImpulse/input.thrust;                        % [s]
input.pcNominal = 70e5;                                                    % [Pa]
input.peNominal = constants.pAmb;                                          % [Pa] (nozzle optimal at sea level)

% Area ratios (wrt At) from which to which there is active cooling
coolARatios = [2, 2];                                            

%% Ballistic characterization
% Vieille's law
[a, aSigma, n, nSigma, R2] = uncertaintyVieille(propellant.ccPressure, propellant.burnRate);

rbNominal = a*input.pcNominal^n;                                           % Burning rate with nominal chamber pressure [mm/s]


%% Ideal thermodynamics & mass sizing
[performanceNom, AtIdeal, AeIdeal] = computeNominalPerformance(input.thrust, propellant.cea.ccTemperature, input.pcNominal, input.peNominal, propellant.cea.gamma, propellant.cea.molarMass, constants);
mPTot = performanceNom.mDot*input.burningTime;                             % Total propellant mass [kg]

%% BATES motor desing
% Grain sizing (using Richard Nakka formulation)
grain = grainDesign(rbNominal, input.burningTime, mPTot, propellant.cea);

% Conical nozzleDesign
nozzle = nozzleDesign(AtIdeal, AeIdeal, grain.dExt0 / 2, ...
                                    "alpha", 15 * pi/180, ...
                                    "beta", 30 * pi/180, "showSummary",true);

%% Combustion chamber / Internal Ballistics (Real Performance)
[t, performance, performanceCEA, grain] = computePerformance(a, n, ...
    propellant, performanceNom, nozzle, grain, constants);

%% Nozzle cooling
% OBTAIN THICKNESS OF TBC BASED ON ENGINEERING CHOICES THAT MAKE FUCKING
% SENSE
design = thermalModelDesign2(propellant, nozzle, performanceCEA, cooling, constants);
out = nozzleThermalModel2(design, input,propellant, nozzle, ...
    performanceCEA, cooling, constants,"showSummary",true);

%% Combustion chamber sizing
cc = combustionChamberDesign(input, performanceNom, grain, nozzle, ...
    propellant, casing, liner, constants);

%% Plots
if plotFlag, plotHollowCylinder(grain); end %#ok<UNRCH>
if plotFlag, plotPerformance(performance, performanceCEA); end %#ok<UNRCH>