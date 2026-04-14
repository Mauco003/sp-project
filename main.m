%% Setup

clear; close all; clc

addpath(genpath("./src"))
plotFlag = 1;

if ~exist(fullfile('.', 'data'), 'dir') || ~exist(fullfile('.', 'data', 'propellant.mat'), 'file')
    data = savePropellantData();
else
    data = load(fullfile('.', 'data', 'propellant.mat'));
end

%% Configuraiton Data

% Generic constants
constants = Constants();

% Requirements / design point
thrust = 100000;                                                            % [N]
totalImpulse = 2.5e6;                                                       % [N s]
burningTime = totalImpulse/thrust;                                          % [s]
pcNominal = 70e5;                                                                  % [Pa]
peNominal = constants.pAmb;                                                        % [Pa] (nozzle optimal at sea level)

% Area ratios (wrt At) from which to which there is active cooling
coolARatios = [2, 2];                                            


%% Ballistic characterization
% Vieille's law
[a, aSigma, n, nSigma, R2] = uncertaintyVieille(data.ccPressure, data.burnRate);

rbNominal = a*pcNominal^n;                                                         % Burning rate with nominal chamber pressure [mm/s]


%% Ideal thermodynamics & mass sizing
[performanceNom, AtIdeal, AeIdeal] = idealPerformance(thrust, data.cea.ccTemperature, pcNominal, peNominal, data.cea.gamma, data.cea.molarMass, constants);

mPTot = performanceNom.mDot*burningTime;                                          % Total propellant mass [kg]


%% BATES motor desing
% Grain sizing (using Richard Nakka formulation)
grain = grainConfiguration(rbNominal, burningTime, mPTot, data.cea);
if plotFlag, plotHollowCylinder(grain); end



% Nozzle geometry (conical) sizing
%
rcc = grain.dExt / 2;                                                       % Chamber radius from grain OD

% Conical nozzleDesign
nozzle = nozzleDesign(AtIdeal, AeIdeal, rcc, ...
                                    "alpha", 15 * pi/180, ...
                                    "beta", 30 * pi/180);


%% Combustion chamber / Internal Ballistics (Real Performance)
[thrust_real, Isp_real, mDot_real] = computePerformance(a, n, ...
    data, performanceNom, nozzle, grain, constants);


