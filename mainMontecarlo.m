% MainMontecarlo - This script runs a Monte Carlo simulation to 
%   analyze the performance of a rocket engine design under 
%   varying conditions. 
%
%   Available parameters:
%      - a: Vieille's law coefficient
%      - n: Vieille's law exponent
%      - O/F: Oxidizer to fuel mass ratio
%      - pa: Ambient pressure
%      - At: Nozzle throat area
%      - Ae: Nozzle exit area
%      - dInt0: Initial grain internal diameter
%      - dExt0: Initial grain external diameter
%      - L0: Initial grain length
%      - mPTot: Total propellant mass
%      - alpha: Nozzle divergent half-angle
%      - beta: Nozzle convergent half-angle



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
pcNominal = 70e5/1.035;                                                                  % [Pa]
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
rcc = grain.dExt0 / 2;                                                       % Chamber radius from grain OD

% Conical nozzleDesign
nozzle = nozzleDesign(AtIdeal, AeIdeal, rcc, ...
                                    "alpha", 15 * pi/180, ...
                                    "beta", 30 * pi/180);


%% Combustion chamber / Internal Ballistics (Real Performance)
[t, performance, grain] = computePerformance(a, n, ...
    data, performanceNom, nozzle, grain, constants);



%% 
% Burning Rate Montecarlo
% histogram of the burning rate considering P nominal and a and n variable
P = pcNominal;

sigma_a_MC = aSigma;
sigma_n_MC = nSigma;

for i = 1:100000
    a_MC = a+sigma_a_MC*randn;
    n_MC = n + sigma_n_MC*randn;
    rb = a_MC*P^n_MC;
    result(i) = rb;
end
histogram(result)

E = mean(result);
E2 = mean(result.^2);
rb_sigma = std(result);
rb_sigma = sqrt(E2-E^2);

%%

% https://jiga.io/articles/aerospace-cnc-machining/
%%
