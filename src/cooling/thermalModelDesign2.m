function design = thermalModelDesign2(propellant, nozzle, performance, cooling, constants, options)

arguments
    propellant struct
    nozzle struct
    performance struct
    cooling struct
    constants Constants
    options.showSummary (1,1) logical = true
end

%% ------------------------------------------------------------------------
% 1) Model configuration
% -------------------------------------------------------------------------
cfg = getModelConfig(nozzle, performance, cooling);

%% ------------------------------------------------------------------------
% 2) Build gas properties at the 5 stations
% -------------------------------------------------------------------------
stations = buildStationsFromCEA(cfg.pc, cfg.stationEps, cfg.stationNames);

%% ------------------------------------------------------------------------
% 3) Search feasible THot_design
% -------------------------------------------------------------------------
best = initializeBest(cfg);

for THot_design = cfg.THot_min : cfg.dTHot : cfg.THot_max

    candidate = evaluateCandidate(THot_design, stations, cfg);

    if candidate.feasible
        design = fillBestWithCandidate(design, candidate, cfg);
        break
    else
        design.lastFail = candidate.failInfo;
    end
end

%% ------------------------------------------------------------------------
% 4) Summary
% -------------------------------------------------------------------------
if options.showSummary
    printSummary(design);
end

end


%% =========================================================================
% Configuration
%% =========================================================================
function cfg = getModelConfig(nozzle, performance, cooling)

% this values are needed for recovery factor and Bartz correlation, as the
% total temperature has been found to increment using the CEA code we hard
% code it as the one in the combustion chamber

cfg.pc = 70e5;                 % [Pa] HARD CODED
cfg.T0 = 2342.92;              % [K]  HARD CODED

cfg.At = nozzle.At;
cfg.Dt = 2*sqrt(cfg.At/pi);

cfg.rCurvature = nozzle.rCurvature;
cfg.cStar = performance.cstar;

% pressure at e=2 * Diam at e = 2 /(2*YS inconel worst value, 600 MPa)* 1.5 SF[m]
inconelYS = 600e6;
cfg.tWall = 65e5*sqrt(2*cfg.At/pi)/inconelYS*1.5;              

cfg.kTBC   = 0.9;                % YSZ [W/m/K]
cfg.kMetal = 24.2;               % INCONEL [W/m/K]

cfg.Thm_max = 1000 + 273.15;   % [K] max of the inconel
cfg.Tcw_max = 150  + 273.15;   % [K] % Boiling point of water

if isfield(cooling, 'Tboil')
    cfg.Tcw_max = cooling.Tboil;
end

cfg.THot_min = 1000;           % [K]
cfg.THot_max = 2400;           % [K]
cfg.dTHot    = 10;             % [K]

cfg.tTBC_min = 100e-6;         % [m]
cfg.tTBC_max = 600e-6;         % [m]

cfg.stationEps   = [2, 1.5, 1, 1.5, 2];
cfg.stationNames = {'Inlet e=2','Conv e=1.5','Throat','Div e=1.5','Exit e=2'};

end


%% =========================================================================
% Build stations from CEA
%% =========================================================================
function stations = buildStationsFromCEA(pc, epsilonList, stationNames)

[~, GasCEA_2]  = getThermoProfileCEA_froz(80, 20, pc/1e5, 2.0, 2.0);
[~, GasCEA_15] = getThermoProfileCEA_froz(80, 20, pc/1e5, 1.5, 1.5);

gammaVals = [GasCEA_2.gamma(1),     GasCEA_15.gamma(1),     GasCEA_2.gamma(2),     GasCEA_15.gamma(3),     GasCEA_2.gamma(3)];
TVals     = [GasCEA_2.T(1),         GasCEA_15.T(1),         GasCEA_2.T(2),         GasCEA_15.T(3),         GasCEA_2.T(3)];
muVals    = [GasCEA_2.viscosity(1), GasCEA_15.viscosity(1), GasCEA_2.viscosity(2), GasCEA_15.viscosity(3), GasCEA_2.viscosity(3)];
cpVals    = [GasCEA_2.cp(1),        GasCEA_15.cp(1),        GasCEA_2.cp(2),        GasCEA_15.cp(3),        GasCEA_2.cp(3)];
prVals    = [GasCEA_2.prandtl(1),   GasCEA_15.prandtl(1),   GasCEA_2.prandtl(2),   GasCEA_15.prandtl(3),   GasCEA_2.prandtl(3)];
machVals  = [GasCEA_2.Mach(1),      GasCEA_15.Mach(1),      GasCEA_2.Mach(2),      GasCEA_15.Mach(3),      GasCEA_2.Mach(3)];
% get also pressure for the thicknes
n = numel(epsilonList);

stations = repmat(struct( ...
    'name', '', ...
    'epsilon', NaN, ...
    'gamma', NaN, ...
    'Tstat', NaN, ...
    'viscosity', NaN, ...
    'cp', NaN, ...
    'prandtl', NaN, ...
    'Mach', NaN), 1, n);

for i = 1:n
    stations(i).name      = stationNames{i};
    stations(i).epsilon   = epsilonList(i);
    stations(i).gamma     = gammaVals(i);
    stations(i).Tstat     = TVals(i);
    stations(i).viscosity = muVals(i);
    stations(i).cp        = cpVals(i);
    stations(i).prandtl   = prVals(i);
    stations(i).Mach      = machVals(i);
end

end


%% =========================================================================
% Evaluate one THot candidate with sequential critical-point logic
%% =========================================================================
function candidate = evaluateCandidate(THot_design, stations, cfg)

nStations = numel(stations);

candidate = initializeCandidate(THot_design);

for i = 1:nStations

    % 1) Size thickness at current station
    [ok, tReq, ~, failReason] = findThicknessAtStation(stations(i), THot_design, cfg);

    if ~ok
        candidate.criticalPhase = i;
        candidate.criticalName  = stations(i).name;
        candidate.failInfo.requiredThickness = tReq;
        candidate.failInfo.station = i;
        candidate.failInfo.reason  = failReason;
        return
    end

    % 2) Check same thickness at next station
    if i < nStations
        nextRes = solveStation(stations(i+1), tReq, THot_design, cfg);

        nextWorks = (nextRes.THot      <= THot_design + 1e-6) && ...
                    (nextRes.THotMetal <= cfg.Thm_max  + 1e-6);

        if nextWorks
            candidate.feasible      = true;
            candidate.tTBC          = tReq;
            candidate.criticalPhase = i;
            candidate.criticalName  = stations(i).name;
            candidate.results       = evaluateAllStations(stations, tReq, THot_design, cfg);
            candidate.failInfo      = emptyFailInfo(THot_design);
            return
        end
    else
        % Last station: if it can be sized, accept it
        candidate.feasible      = true;
        candidate.tTBC          = tReq;
        candidate.criticalPhase = i;
        candidate.criticalName  = stations(i).name;
        candidate.results       = evaluateAllStations(stations, tReq, THot_design, cfg);
        candidate.failInfo      = emptyFailInfo(THot_design);
        return
    end
end

candidate.failInfo.reason = 'No critical station found with sequential logic';

end


%% =========================================================================
% Find thickness such that THot = THot_target at one station
% Cold-side wall temperature is fixed to Tcw_max
%% =========================================================================
function [ok, tReq, resReq, failReason] = findThicknessAtStation(station, THot_target, cfg)

ok = false;
tReq = NaN;
resReq = struct([]);
failReason = '';

resMin = solveStation(station, cfg.tTBC_min, THot_target, cfg);
resMax = solveStation(station, cfg.tTBC_max, THot_target, cfg);

fMin = resMin.THot - THot_target;
fMax = resMax.THot - THot_target;

% THot increases with thickness for fixed cold-side temperature
if fMin > 0
    tReq = cfg.tTBC_min;
    resReq = resMin;
    failReason = 'THot_target below minimum achievable THot at tMin';
    return
end

if fMax < 0
    tReq = cfg.tTBC_max;
    resReq = resMax;
    failReason = 'THot_target requires thickness larger than tTBC_max_search';
    return
end

a = cfg.tTBC_min;
b = cfg.tTBC_max;
resB = resMax;

for iter = 1:50
    m = 0.5*(a+b);
    resM = solveStation(station, m, THot_target, cfg);
    fM = resM.THot - THot_target;

    if fM < 0
        a = m;
    else
        b = m;
        resB = resM;
    end
end

tReq = b;
resReq = resB;

if resReq.THotMetal > cfg.Thm_max + 1e-6
    failReason = 'THot_target matched but THotMetal exceeds Thm_max';
    return
end

ok = true;

end


%% =========================================================================
% Recompute all stations with the chosen thickness
%% =========================================================================
function results = evaluateAllStations(stations, tChosen, THot_ref, cfg)

nStations = numel(stations);

results = repmat(struct( ...
    'name', '', ...
    'epsilon', NaN, ...
    'gamma', NaN, ...
    'mach', NaN, ...
    'Tstat', NaN, ...
    'Taw', NaN, ...
    'hGas', NaN, ...
    'qDot', NaN, ...
    'THot', NaN, ...
    'THotMetal', NaN, ...
    'TColdMetal', NaN, ...
    'tRequired', tChosen), 1, nStations);

for i = 1:nStations
    res = solveStation(stations(i), tChosen, THot_ref, cfg);

    results(i).name       = stations(i).name;
    results(i).epsilon    = stations(i).epsilon;
    results(i).gamma      = stations(i).gamma;
    results(i).mach       = stations(i).Mach;
    results(i).Tstat      = stations(i).Tstat;
    results(i).Taw        = res.Taw;
    results(i).hGas       = res.hGas;
    results(i).qDot       = res.qDot;
    results(i).THot       = res.THot;
    results(i).THotMetal  = res.THotMetal;
    results(i).TColdMetal = res.TColdMetal;
    results(i).tRequired  = tChosen;
end

end


%% =========================================================================
% Solve one station for a fixed thickness
% Boundary condition: TColdMetal = Tcw_max
%% =========================================================================
function res = solveStation(station, tTBC, THot_ref, cfg)

gamma = station.gamma;
Tstat = station.Tstat;
muGas = station.viscosity;
cpGas = station.cp;
PrGas = station.prandtl;
mach  = station.Mach;
epsi  = station.epsilon;

r = recoveryFactor(gamma, mach, PrGas);
Taw = cfg.T0 * r;

Rcond = tTBC/cfg.kTBC + cfg.tWall/cfg.kMetal;

THot_guess = min(max(THot_ref, cfg.Tcw_max + 1e-6), Taw - 1e-6);

for iter = 1:50
    hGas = bartzCorrelation(cfg.pc, cfg.cStar(1), cfg.Dt, cfg.rCurvature, epsi, ...
                            muGas, cpGas, PrGas, THot_guess, cfg.T0, Tstat, 0.6);

    qDot = (Taw - cfg.Tcw_max) / (1/hGas + Rcond);

    THot_new      = cfg.Tcw_max + qDot*Rcond;
    THotMetal_new = cfg.Tcw_max + qDot*cfg.tWall/cfg.kMetal;

    if abs(THot_new - THot_guess) < 1e-8
        THot_guess = THot_new;
        break
    end

    THot_guess = 0.5*THot_guess + 0.5*THot_new;
end

res.Taw        = Taw;
res.hGas       = hGas;
res.qDot       = qDot;
res.THot       = THot_new;
res.THotMetal  = THotMetal_new;
res.TColdMetal = cfg.Tcw_max;

end


%% =========================================================================
% Best struct helpers
%% =========================================================================
function best = initializeBest(cfg)

best.found = false;
best.results = struct([]);

best.THot = NaN;
best.tTBC = NaN;
best.tWall = cfg.tWall;

best.kTBC = cfg.kTBC;
best.kMetal = cfg.kMetal;

best.Thm_max = cfg.Thm_max;
best.Tcw_max = cfg.Tcw_max;

best.criticalPhase = NaN;
best.criticalName = '';

best.maxTHot = NaN;
best.minTHot = NaN;
best.maxTHotMetal = NaN;
best.minTHotMetal = NaN;
best.maxTColdMetal = NaN;
best.minTColdMetal = NaN;
best.maxQDot = NaN;

best.qDotCritical = NaN;
best.THotCritical = NaN;
best.THotMetalCritical = NaN;
best.TColdMetalCritical = NaN;

best.lastFail = emptyFailInfo(NaN);

end


function best = fillBestWithCandidate(best, candidate, cfg)

results = candidate.results;

best.found = true;
best.results = results;

best.THot = candidate.failInfo.THot;
best.tTBC = candidate.tTBC;

best.criticalPhase = candidate.criticalPhase;
best.criticalName  = candidate.criticalName;

allTH  = [results.THot];
allTHM = [results.THotMetal];
allTCW = [results.TColdMetal];
allQ   = [results.qDot];

best.maxTHot       = maxNoNaN(allTH);
best.minTHot       = minNoNaN(allTH);
best.maxTHotMetal  = maxNoNaN(allTHM);
best.minTHotMetal  = minNoNaN(allTHM);
best.maxTColdMetal = maxNoNaN(allTCW);
best.minTColdMetal = minNoNaN(allTCW);
best.maxQDot       = maxNoNaN(allQ);

iCrit = candidate.criticalPhase;
best.qDotCritical       = results(iCrit).qDot;
best.THotCritical       = results(iCrit).THot;
best.THotMetalCritical  = results(iCrit).THotMetal;
best.TColdMetalCritical = results(iCrit).TColdMetal;

best.lastFail = emptyFailInfo(NaN);
best.lastFail.THot = candidate.failInfo.THot;

best.tWall   = cfg.tWall;
best.kTBC    = cfg.kTBC;
best.kMetal  = cfg.kMetal;
best.Thm_max = cfg.Thm_max;
best.Tcw_max = cfg.Tcw_max;

end


%% =========================================================================
% Candidate helpers
%% =========================================================================
function candidate = initializeCandidate(THot_design)

candidate.feasible = false;
candidate.tTBC = NaN;
candidate.criticalPhase = NaN;
candidate.criticalName = '';
candidate.results = struct([]);
candidate.failInfo = emptyFailInfo(THot_design);

end


function failInfo = emptyFailInfo(THot_design)

failInfo = struct( ...
    'THot', THot_design, ...
    'requiredThickness', NaN, ...
    'station', NaN, ...
    'reason', '');

end


%% =========================================================================
% Printing
%% =========================================================================
function printSummary(best)

fprintf('\n');
fprintf('====================================================\n');
fprintf('         THERMAL MODEL - JACKET DESIGN              \n');
fprintf('====================================================\n');
fprintf('Found feasible solution? : %d\n', best.found);

if ~best.found
    fprintf('No feasible solution found in the explored range.\n');
    fprintf('Last checked THot        : %.3f K\n', best.lastFail.THot);
    fprintf('Required thickness       : %.6e m\n', best.lastFail.requiredThickness);
    fprintf('Failure station          : %d\n', best.lastFail.station);
    fprintf('Failure reason           : %s\n', best.lastFail.reason);
    fprintf('====================================================\n\n');
    return
end

fprintf('THot design              : %.3f K\n', best.THot);
fprintf('tTBC                     : %.6e m\n', best.tTBC);
fprintf('tWall                    : %.6e m\n', best.tWall);
fprintf('Critical phase           : %d (%s)\n', best.criticalPhase, best.criticalName);
fprintf('Critical qDot            : %.3f W/m^2\n', best.qDotCritical);
fprintf('Critical THot            : %.3f K\n', best.THotCritical);
fprintf('Critical THotMetal       : %.3f K\n', best.THotMetalCritical);
fprintf('Critical TColdMetal      : %.3f K\n', best.TColdMetalCritical);
fprintf('THot local min / max     : %.3f / %.3f K\n', best.minTHot, best.maxTHot);
fprintf('THotMetal min / max      : %.3f / %.3f K\n', best.minTHotMetal, best.maxTHotMetal);
fprintf('TColdMetal min / max     : %.3f / %.3f K\n', best.minTColdMetal, best.maxTColdMetal);
fprintf('Max qDot                 : %.3f W/m^2\n', best.maxQDot);
fprintf('====================================================\n\n');

end


%% =========================================================================
% Utilities
%% =========================================================================
function val = maxNoNaN(x)
x = x(~isnan(x));
if isempty(x)
    val = NaN;
else
    val = max(x);
end
end

function val = minNoNaN(x)
x = x(~isnan(x));
if isempty(x)
    val = NaN;
else
    val = min(x);
end
end