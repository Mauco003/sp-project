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
hoopStress = 120e6;       % update these w real vals
safetyFactor = 1.5;     % update w real values
pc = 70;

casing = struct();
casing.thickness = (pc*d_grain)/(2*hoopStress) * safetyFactor;
casing.thermalConductivity = 1;
casing.TMax = 100;
casing.cost = 1;        % [cost per kg]
casing.density = 1;

%% LINER
liner = struct();
liner.thickness = 1;
liner.thermalConductivity = 1;
liner.regressionRate = 1;
liner.cost = 1;         % [cost per kg]
liner.density = 1;

%% OPTIONS
t_low = 1e-3;
t_high = 10e-3;
iter = 100;
wMass = 1;
wCost = 1;

%% ITERATIVE FOR LOOP FOR THICKNESSES
% for i = 1:iter
% 
%     % guess for liner thickness
%     t_mid = (t_low + t_high)/2;
%     liner.thickness = t_mid;
% 
%     % check if the thickness is optimal
%     results = combustionChamber(casing, liner, nozzle, propellant, performance, constants);
% 
%     % adjust value of liner thickness
%     if results.casingSurvives
%         t_high = t_mid;
%     else
%         t_low = t_mid;
%     end
% end
% 
% t_opt = t_high;



tL_vals = linspace(1e-3, 20e-3, 40);   % liner thickness [m]
tC_vals = linspace(1e-3, 20e-3, 40);   % casing thickness [m]

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