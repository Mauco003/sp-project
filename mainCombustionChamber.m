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

% geometry of grain
d_grain = 0.89;
r_grain = d_grain/2;
l_grain = 1.58;

% Area ratios (wrt At) from which to which there is active cooling
coolARatios = [2, 2];                                            

%% Ballistic characterization
% Vieille's law
[a, aSigma, n, nSigma, R2] = uncertaintyVieille(propellant.ccPressure, propellant.burnRate);

rbNominal = a*input.pcNominal^n;                                                         % Burning rate with nominal chamber pressure [mm/s]


%% Ideal thermodynamics & mass sizing
[performanceNom, AtIdeal, AeIdeal] = computeNominalPerformance(input.thrust, propellant.cea.ccTemperature, input.pcNominal, input.peNominal, propellant.cea.gamma, propellant.cea.molarMass, constants);
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
hoopStress = 460e6;       % update these w real vals
safetyFactor = 1.5;     % update w real values
pc = 70e5;

casing = struct();
casing.thickness = (pc*d_grain)/(2*hoopStress) * safetyFactor;
casing.thermalConductivity = 42.7;
casing.TMax = 1432;
casing.cost = 0.8;        % [cost per kg]
casing.density = 7850;   % kg/m^3

%% LINER
liner = struct();
liner.thickness = 5e-3;
liner.thermalConductivity = 0.125;
liner.regressionRate = 0.225;         % kg/s*m^2
liner.cost = 7;         % [cost per kg]
liner.density = 1208;   % kg/m^3

%% OPTIONS
t_low = 1e-3;
t_high = 10e-3;
iter = 100;
wMass = 1;
wCost = 1;

%% optimization loop
% ranges to search
tL_vals = linspace(1e-5, 5e-2, 100);   % liner thickness [m]
tC_vals = linspace(casing.thickness, 20e-2, 100);   % casing thickness [m]

best.J = 1000000000000;

for i = 1:length(tL_vals)
    for j = 1:length(tC_vals)

        % Set thicknesses
        liner.thickness  = tL_vals(i);
        casing.thickness = tC_vals(j);

        % Run thermal model
        results = combustionChamber(casing, liner, nozzle, propellant, performance, constants, ...
            "makePlot", false, "showSummary", false);

        % Check survival constraint
        if results.casingSurvives && results.linerSurvives

            % Compute cost + mass
            [J, massTotal, costTotal] = ccCost(liner, casing, wMass, wCost);

            % Store best solution
            if J < best.J
                best.J = J;
                best.t_liner = liner.thickness;
                best.t_casing = casing.thickness;
                best.mass = massTotal;
                best.cost = costTotal;
                best.results = results;
            end

        end
    end
end



% liner nozzle thickness
t_burn = 25; 
t_liner_nozzle = (liner.regressionRate / liner.density * t_burn)*2;