%% Setup
clear; close all; clc

addpath(genpath("."))
plotFlag = 0;

%% Input Data

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

% Empirical data
% Pressure [bar]
pressureData = [10.1; 9.7; 10.3; 30.2; 31.0; 29.8; 50.2; 51.0; ...
            50.3; 70.2; 69.0; 69.8; 91.2; 89.1; 89.0];

% Burning rate [mm/s]
burningRateData = [4.0; 3.8; 4.1; 5.6; 6.0; 5.7; 7.0; 7.2; ...
                7.1; 8.4; 8.3; 8.6; 8.8; 9.0; 9.2];

% from NASA CEA output 80% AP - 20% HTPB 
gamma = 1.2386*1.00032;                                                     % Specific heat ratio [-]
Tc = 2342.92;                                                               % Chamber temperature [K]
molarMass = 22.087;

%% 

% Vieille's law
[a, aSigma, n, nSigma, R2] = Uncertainty(pressureData, burningRateData);


% Conical nozzleDesign
[nozzle, propPerf] = nozzleDesign(thrust, Tc, pc, pe, gamma, molarMass, constants, ...
                                    "alpha", 15 * pi/180, ...
                                    "beta", 30 * pi/180);

mPTot = propPerf.mDot*burningTime;

% Combustion chamber
