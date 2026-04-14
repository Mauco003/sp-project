%% Setup

clear; close all; clc

addpath(genpath("./src"))
plotFlag = 1;


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

% from NASA CEA output- composition: 80% AP - 20% HTPB
gamma = 1.2386*1.00032;                                                     % Specific heat ratio (gammas * (dLV/dLP)t) [-]
Tc = 2342.92;                                                               % Chamber temperature [K]
rhoAp = 1950;                                                               % Ammonium perchlorate density [kg/m^3]
rhoHTPB = 913;                                                              % HTPB density [kg/m^3]
molarMassProp = 22.087;
rhoP = 1/(0.8/rhoAp +  0.2/rhoHTPB);                                        % Propellant density [kg/m^3]


%% 

% Vieille's law
[a, aSigma, n, nSigma, R2] = Uncertainty(pressureData, burningRateData);


% Conical nozzleDesign
[nozzle, propPerf] = nozzleDesign(thrust, Tc, pc, pe, gamma, molarMassProp, ...
                                    constants, ...
                                    "alpha", 15 * pi/180, ...
                                    "beta", 30 * pi/180);

mPTot = propPerf.mDot*burningTime;                                          % Total propellant mass [kg]

rb = a*(pc*1e-5)^n;                                                         % Burning rate with nominal chamber pressure [mm/s]
web = rb*burningTime;
Vprop = mPTot/rhoP;

% BATES motor design (using Richard Nakka formulation)
dInt = internalDiameter(Vprop, web*1e-3); % m
[dExt, L0] = grainConfiguration(dInt, web*1e-3); % m

if plotFlag, plotHollowCylinder(dExt, dInt, L0); end

% Combustion chamber
% [t, p, rb] = computeBurn(a, n, rhoP*1e-6, propPerf.cstar, dExt, dInt, L0, nozzle.At*1e6);

