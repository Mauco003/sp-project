function best = thermalModelMaggi(propellant, nozzle, performance, cooling, constants, options)

arguments
    propellant struct
    nozzle struct
    performance struct
    cooling struct
    constants Constants
    options.showSummary (1,1) logical = true
end

%% ------------------------------------------------------------------------
% 1) Inputs
% -------------------------------------------------------------------------
pc = 70e5;                 % [Pa]
T0 = 2342.92;              % HARD CODED DONT CHANGE PLEASE

At = nozzle.At;
Dt = 2*sqrt(At/pi);

rCurvature = nozzle.rCurvature;   % [m]
cStar = performance.cstar;

tWall = 1e-3;              % [m]

kTBC   = 1;                % [W/m/K]
kMetal = 22;               % [W/m/K]

Thm_max = 1800 + 273.15;   % [K]
Tcw_max = 180 + 273.15;    % [K]

epsilon = [2, 1.5, 1, 1.5, 2];
stationNames = {'Inlet e=2','Conv e=1.5','Throat','Div e=1.5','Exit e=2'};

%% ------------------------------------------------------------------------
% 2) Build 5 stations from CEA
%    Assumption:
%    getThermoProfileCEA(..., eConv, eDiv) returns 3 stations:
%    [station at eConv, throat, station at eDiv]
% -------------------------------------------------------------------------
[~, GasCEA_2]   = getThermoProfileCEA(80, 20, pc/1e5, 2.0, 2.0);
[~, GasCEA_15]  = getThermoProfileCEA(80, 20, pc/1e5, 1.5, 1.5);

% Build 5-point station data
Gas5.gamma     = [GasCEA_2.gamma(1),     GasCEA_15.gamma(1),     GasCEA_2.gamma(2),     GasCEA_15.gamma(3),     GasCEA_2.gamma(3)];
Gas5.T         = [GasCEA_2.T(1),         GasCEA_15.T(1),         GasCEA_2.T(2),         GasCEA_15.T(3),         GasCEA_2.T(3)];
Gas5.viscosity = [GasCEA_2.viscosity(1), GasCEA_15.viscosity(1), GasCEA_2.viscosity(2), GasCEA_15.viscosity(3), GasCEA_2.viscosity(3)];
Gas5.cp        = [GasCEA_2.cp(1),        GasCEA_15.cp(1),        GasCEA_2.cp(2),        GasCEA_15.cp(3),        GasCEA_2.cp(3)];
Gas5.prandtl   = [GasCEA_2.prandtl(1),   GasCEA_15.prandtl(1),   GasCEA_2.prandtl(2),   GasCEA_15.prandtl(3),   GasCEA_2.prandtl(3)];
Gas5.Mach      = [GasCEA_2.Mach(1),      GasCEA_15.Mach(1),      GasCEA_2.Mach(2),      GasCEA_15.Mach(3),      GasCEA_2.Mach(3)];

%% ------------------------------------------------------------------------
% 3) Search on THot and tTBC
% -------------------------------------------------------------------------
found = false;
results = struct([]);

for THot = 1500 : 10 : 2400
    for tTBC = 100e-6 : 50e-6 : 600e-6

        [stationOK, resultsCandidate] = evaluateCandidate( ...
            Gas5, epsilon, stationNames, ...
            pc, cStar, Dt, rCurvature, ...
            T0, THot, tTBC, tWall, ...
            kTBC, kMetal, Thm_max, Tcw_max, cooling);

        if stationOK
            found = true;
            results = resultsCandidate;
            break
        end
    end

    if found
        break
    end
end

%% ------------------------------------------------------------------------
% 4) If no feasible solution found, keep last evaluated candidate
% -------------------------------------------------------------------------
if ~found
    results = resultsCandidate;
end

%% ------------------------------------------------------------------------
% 5) Outputs
% -------------------------------------------------------------------------
best.results = results;

best.THot = THot;
best.tTBC = tTBC;
best.tWall = tWall;

best.kTBC = kTBC;
best.kMetal = kMetal;

best.Thm_max = Thm_max;
best.Tcw_max = Tcw_max;

best.found = found;

allTHM = [results.THotMetal];
allTCW = [results.TColdMetal];
allQ   = [results.qDot];

best.maxTHotMetal = max(allTHM(:));
best.maxTColdMetal = max(allTCW(:));
best.maxQDot = max(allQ(:));

%% ------------------------------------------------------------------------
% 6) Summary
% -------------------------------------------------------------------------
if options.showSummary
    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('              THERMAL MODEL MAGGI                   \n');
    fprintf('====================================================\n');
    fprintf('Found feasible solution? : %d\n', best.found);
    fprintf('THot                     : %.3f K\n', best.THot);
    fprintf('tTBC                     : %.6e m\n', best.tTBC);
    fprintf('tWall                    : %.6e m\n', best.tWall);
    fprintf('Max THotMetal            : %.3f K\n', best.maxTHotMetal);
    fprintf('Max TColdMetal           : %.3f K\n', best.maxTColdMetal);
    fprintf('Max qDot                 : %.3f W/m^2\n', best.maxQDot);
    fprintf('====================================================\n');

    for i = 1:numel(results)
        fprintf('\n');
        fprintf('--- %s ---\n', results(i).name);
        fprintf('epsilon                  : %.3f\n', results(i).epsilon);
        fprintf('gamma                    : %.4f\n', results(i).gamma);
        fprintf('Mach                     : %.4f\n', results(i).mach);
        fprintf('Tstat                    : %.3f K\n', results(i).Tstat);
        fprintf('Taw                      : %.3f K\n', results(i).Taw);
        fprintf('hGas                     : %.3f W/m^2/K\n', results(i).hGas);
        fprintf('qDot                     : %.3f W/m^2\n', results(i).qDot);
        fprintf('THotMetal                : %.3f K\n', results(i).THotMetal);
        fprintf('TColdMetal               : %.3f K\n', results(i).TColdMetal);
    end

    fprintf('\n');
end

end


%% ========================================================================
% Local evaluator for one candidate (THot, tTBC)
%% ========================================================================
function [stationOK, results] = evaluateCandidate( ...
    Gas5, epsilon, stationNames, ...
    pc, cStar, Dt, rCurvature, ...
    T0, THot, tTBC, tWall, ...
    kTBC, kMetal, Thm_max, Tcw_max, cooling)

nStations = numel(epsilon);
results = repmat(struct( ...
    'name', '', ...
    'epsilon', NaN, ...
    'gamma', NaN, ...
    'mach', NaN, ...
    'Tstat', NaN, ...
    'Taw', NaN, ...
    'hGas', NaN, ...
    'qDot', NaN, ...
    'THotMetal', NaN, ...
    'TColdMetal', NaN), 1, nStations);

stationOK = true;

for i = 1:nStations
    gamma = Gas5.gamma(i);
    Tstat = Gas5.T(i);
    muGas = Gas5.viscosity(i);
    cpGas = Gas5.cp(i);
    PrGas = Gas5.prandtl(i);
    mach  = Gas5.Mach(i);

    hGas = bartzCorrelation(pc, cStar, Dt, rCurvature, epsilon(i), muGas, cpGas, PrGas, THot, T0, Tstat, 0.6);

    r = recoveryFactor(gamma, mach, PrGas);
    Taw = T0 * r;

    qDot = hGas * (Taw - THot);

    THotMetal = THot - qDot*tTBC/kTBC;
    TColdMetal = THotMetal - qDot*tWall/kMetal;

    results(i).name = stationNames{i};
    results(i).epsilon = epsilon(i);
    results(i).gamma = gamma;
    results(i).mach = mach;
    results(i).Tstat = Tstat;
    results(i).Taw = Taw;
    results(i).hGas = hGas;
    results(i).qDot = qDot;
    results(i).THotMetal = THotMetal;
    results(i).TColdMetal = TColdMetal;

    if any(THotMetal > Thm_max) || any(TColdMetal > Tcw_max)
        stationOK = false;
        return
    end
end

end