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
pc = 70e5;                                                                  % [Pa]
pe = constants.pAmb;                                                        % [Pa] (nozzle optimal at sea level)

% Area ratios (wrt At) from which to which there is active cooling
coolARatios = [2, 2];                                            

% Vieille's law
[a, aSigma, n, nSigma, R2] = uncertaintyVieille(data.ccPressure, data.burnRate);


% Conical nozzleDesign
[nozzle, propPerf] = nozzleDesign(thrust, data.cea.ccTemperature, pc, pe, data.cea.gamma, data.cea.molarMass, ...
                                    constants, ...
                                    "alpha", 15 * pi/180, ...
                                    "beta", 30 * pi/180);

mPTot = propPerf.mDot*burningTime;                                          % Total propellant mass [kg]

rb = a*pc^n;                                                         % Burning rate with nominal chamber pressure [mm/s]
web = rb*burningTime;
Vprop = mPTot/data.cea.rhoP;

% BATES motor design (using Richard Nakka formulation)
% dInt = internalDiameter(Vprop, web*1e-3); % m
[dExt, L0] = grainConfiguration(Vprop, web); % m

if plotFlag, plotHollowCylinder(dExt, dInt, L0); end

% Combustion chamber
% [t, p, rb] = computeBurn(a, n, data.cea.rhoP*1e-6, propPerf.cstar, dExt, dInt, L0, nozzle.At*1e6);

