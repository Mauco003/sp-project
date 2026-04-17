%% MAIN TO RUN COMBUSTION CHAMBER FUNCTION -------------------------------
% xx
%
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
[performanceNom, AtIdeal, AeIdeal] = performanceNomCalc(input.thrust, propellant.cea.ccTemperature, input.pcNominal, input.peNominal, propellant.cea.gamma, propellant.cea.molarMass, constants);
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


%% CASING
casing = struct();
casing.thickness = 1;
casing.thermalConductivity = 1;
casing.TMax = 100;

%% LINER
liner = struct();
liner.thickness = 1;
liner.thermalConductivity = 1;
liner.regressionRate = 1;

%% OPTIONS
t_low = 1e-3;
t_high = 10e-3;


%% ITERATIVE FOR LOOP FOR THICKNESSES
for i = 1:30
    t_mid = (t_low + t_high)/2;
    liner.thickness = t_mid;

    results = combustionChamber(casing, liner, nozzle, propellant, performance, constants);

    if results.casingSurvives
        t_high = t_mid;
    else
        t_low = t_mid;
    end
end

t_opt = t_high;