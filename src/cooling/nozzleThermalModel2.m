function out = nozzleThermalModel2(design, input, propellant, nozzle, performance, cooling, constants, options)

arguments
    design
    input
    propellant                struct
    nozzle                    struct
    performance               struct
    cooling                   struct
    constants                 Constants
    options.n           (1,1) double  = 10
    options.areaLimits  (2,1) double  = [2; 2]
    options.makePlot    (1,1) logical = true
    options.showSummary (1,1) logical = true
    options.mDotList    (1,:) double  = linspace(1,10,10)
end

%% ------------------------------------------------------------------------
% 0) Fixed design values from design model
% -------------------------------------------------------------------------
tWall   = design.tWall;
tTBC    = design.tTBC*1.2;
kWall   = design.kMetal;
kTBC    = design.kTBC;
TwallMaxCoolantSide = design.Tcw_max;   % usually 150 C in K

%% ------------------------------------------------------------------------
% 1) Fixed coolant data
% -------------------------------------------------------------------------
cpWater  = 4180;
rhoWater = 997;
muWater  = 0.89e-3;
kWater   = 0.60;
PrWater  = cpWater * muWater / kWater;

waterT0       = 18 + 273.15;
waterPressure = 10e5;

mDotList = options.mDotList(:).';
nCases   = numel(mDotList);

%% ------------------------------------------------------------------------
% 2) Gas/chamber data
% -------------------------------------------------------------------------
R = constants.R / propellant.cea.molarMass;

gamma = propellant.cea.gamma;
pc    = propellant.ccPressure;
cStar = performance.cstar;
muGas = propellant.cea.mu;
cpGas = R * gamma/(gamma-1);
PrGas = muGas*cpGas / propellant.cea.k;

%% ------------------------------------------------------------------------
% 3) Nozzle geometry and mesh
% -------------------------------------------------------------------------
At    = nozzle.At;
alpha = nozzle.alpha;
beta  = nozzle.beta;

n = options.n;
N = 2*n - 1;
areaRatioCool = options.areaLimits;

rt = sqrt(At/pi);
rLim = sqrt(areaRatioCool*At/pi);

Lconv = (rLim(1) - rt)/tan(beta);
Ldiv  = (rLim(2) - rt)/tan(alpha);

xConv = linspace(-Lconv, 0, n)';
xDiv  = linspace(0, Ldiv, n)';

x = [xConv; xDiv(2:end)];
r = zeros(size(x));

r(x <= 0) = rt - x(x <= 0)*tan(beta);
r(x > 0)  = rt + x(x > 0)*tan(alpha);

A = pi*r.^2;
epsilon = A/At;

dA = 2*pi*(sqrt(diff(x).^2 + diff(r).^2)) .* 0.5 .* (r(1:end-1) + r(2:end));

%% ------------------------------------------------------------------------
% 4) CEA interpolation on mesh
% -------------------------------------------------------------------------
[machCEA, TstatCEA, gammaCEA, cpCEA, muCEA, PrCEA] = ...
    buildCEAProfileOnMeshSimple(epsilon, n, input.pcNominal, muGas, PrGas);

%% ------------------------------------------------------------------------
% 5) Run all requested mass flows
% -------------------------------------------------------------------------
T_CEA_all        = cell(1,nCases);
q_CEA_all        = cell(1,nCases);
hGas_CEA_all     = cell(1,nCases);
hWater_CEA_all   = cell(1,nCases);
T0_CEA_all       = cell(1,nCases);
maxQCEAList      = zeros(1,nCases);

TwaterOutList    = NaN(1,nCases);
maxTcoldList     = NaN(1,nCases);
validMdotMask    = false(1,nCases);

for j = 1:nCases

    mDot = mDotList(j);

    T_CEA      = zeros(N,5);   % [Taw, Thot, T_TBC_metal, Tcold, Twater]
    q_CEA      = zeros(N,1);
    hGas_CEA   = zeros(N,1);
    hWater_CEA = zeros(N,1);
    T0_CEA     = zeros(N,1);

    T_CEA(1,5) = waterT0;

    for i = 1:N
        T0_CEA(i)   = TstatCEA(i) * (1 + 0.5*(gammaCEA(i)-1)*machCEA(i)^2);
        T_CEA(i,1)  = T0_CEA(i) * recoveryFactor(gammaCEA(i), machCEA(i), PrCEA(i));

        [hWater_i, ~] = coolantHTCfromMassFlow(mDot, rhoWater, muWater, kWater, PrWater);
        hWater_CEA(i) = hWater_i;

        [q_i, hGas_i, Thot_i, Tmid_i, Tcold_i] = ...
            solveStationIterative( ...
                T_CEA(i,1), T_CEA(i,5), ...
                input.pcNominal, performance.cstar, 2*rt, nozzle.rCurvature, epsilon(i), ...
                muCEA(i), cpCEA(i), PrCEA(i), ...
                T0_CEA(i), TstatCEA(i), ...
                hWater_i, tWall, tTBC, kWall, kTBC);

        q_CEA(i)    = q_i;
        hGas_CEA(i) = hGas_i;
        T_CEA(i,2)  = Thot_i;
        T_CEA(i,3)  = Tmid_i;
        T_CEA(i,4)  = Tcold_i;

        if i < N
            T_CEA(i+1,5) = T_CEA(i,5) + q_CEA(i)*dA(i)/(cpWater*mDot);
        end
    end

    T_CEA_all{j}      = T_CEA;
    q_CEA_all{j}      = q_CEA;
    hGas_CEA_all{j}   = hGas_CEA;
    hWater_CEA_all{j} = hWater_CEA;
    T0_CEA_all{j}     = T0_CEA;

    maxQCEAList(j)   = max(q_CEA);
    maxTcoldList(j)  = max(T_CEA(:,4));

    % Validity criterion:
    % the maximum coolant-side wall temperature must not exceed Tboil limit
    if maxTcoldList(j) <= TwallMaxCoolantSide + 1e-9
        validMdotMask(j) = true;
        TwaterOutList(j) = T_CEA(end,5);
    else
        validMdotMask(j) = false;
        TwaterOutList(j) = NaN;   % invalid so plots do not break
    end
end

%% ------------------------------------------------------------------------
% 6) Select reference case for detailed plots
% -------------------------------------------------------------------------
validIdx = find(validMdotMask);

if isempty(validIdx)
    idxRef = [];
    warning('No valid mass flow found: all cases exceed the coolant-side wall temperature limit.');
else
    [~, kRef] = min(abs(mDotList(validIdx) - 3));
    idxRef = validIdx(kRef);
end

if ~isempty(idxRef)
    mDotRef        = mDotList(idxRef);
    T_CEA_ref      = T_CEA_all{idxRef};
    q_CEA_ref      = q_CEA_all{idxRef};
    hGas_CEA_ref   = hGas_CEA_all{idxRef};
    hWater_CEA_ref = hWater_CEA_all{idxRef};
    T0_CEA_ref     = T0_CEA_all{idxRef};
else
    mDotRef        = NaN;
    T_CEA_ref      = [];
    q_CEA_ref      = [];
    hGas_CEA_ref   = [];
    hWater_CEA_ref = [];
    T0_CEA_ref     = [];
end

%% ------------------------------------------------------------------------
% 7) Outputs
% -------------------------------------------------------------------------
out = struct();

out.x = x;
out.r = r;
out.A = A;
out.epsilon = epsilon;
out.dA = dA;

out.tWall = tWall;
out.tTBC  = tTBC;

out.mDotList = mDotList;
out.mDotRef  = mDotRef;

out.waterOutletTemperatureList = TwaterOutList;
out.maxQCEAList = maxQCEAList;
out.maxTcoldList = maxTcoldList;
out.validMdotMask = validMdotMask;
out.coolantSideWallTemperatureLimit = TwallMaxCoolantSide;

out.TCEA = T_CEA_ref;
out.qCEA = q_CEA_ref;
out.hGasCEA = hGas_CEA_ref;
out.hWaterCEA = hWater_CEA_ref;
out.T0_CEA = T0_CEA_ref;

out.TCEA_all = T_CEA_all;
out.qCEA_all = q_CEA_all;
out.hGasCEA_all = hGas_CEA_all;

out.machCEA   = machCEA;
out.TstatCEA  = TstatCEA;
out.gammaCEA  = gammaCEA;
out.cpCEA     = cpCEA;

out.design = design;
out.waterInletTemperature = waterT0;
out.waterPressure = waterPressure;

%% ------------------------------------------------------------------------
% 8) Summary
% -------------------------------------------------------------------------
if options.showSummary
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('            NOZZLE THERMAL MODEL - DISCRETIZED              \n');
    fprintf('============================================================\n');
    fprintf('tWall                         : %.6e m\n', tWall);
    fprintf('tTBC                          : %.6e m\n', tTBC);
    fprintf('Water inlet temp              : %.2f C\n', waterT0 - 273.15);
    fprintf('Coolant-side wall temp limit  : %.2f C\n', TwallMaxCoolantSide - 273.15);

    if ~isempty(idxRef)
        fprintf('Reference mDot                : %.4f kg/s\n', mDotRef);
    else
        fprintf('Reference mDot                : none (no valid case)\n');
    end

    fprintf('\n');
    fprintf('mDot [kg/s]    Tout [C]    Max q [W/m^2]    Max Tcold [C]    Valid\n');
    for j = 1:nCases
        if isnan(TwaterOutList(j))
            toutStr = '   NaN  ';
        else
            toutStr = sprintf('%8.2f', TwaterOutList(j) - 273.15);
        end
        fprintf('%10.4f    %s    %14.4e    %12.2f    %d\n', ...
            mDotList(j), toutStr, maxQCEAList(j), maxTcoldList(j) - 273.15, validMdotMask(j));
    end
    fprintf('============================================================\n\n');
end

%% ------------------------------------------------------------------------
% 9) Plots
%% -------------------------------------------------------------------------
if options.makePlot

    saveDir = 'C:\Users\celia\Desktop\figuras_cooling';
    legendFontSize = 18;

    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    % 1) Outlet water temperature vs mass flow
    fig = figure('Color', 'w', 'Position', [100, 100, 850, 520], ...
                 'Name', 'water-outlet-temperature-vs-mass-flow');
    hold on;

    validIdxPlot   = find(validMdotMask);
    invalidIdxPlot = find(~validMdotMask);

    if ~isempty(validIdxPlot)
        plot(mDotList(validIdxPlot), TwaterOutList(validIdxPlot) - 273.15, '-', ...
            'Color', [0.65 0.65 0.65], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot(mDotList(validIdxPlot), TwaterOutList(validIdxPlot) - 273.15, 'o', ...
            'MarkerEdgeColor', [0.2 0.2 0.2], ...
            'MarkerFaceColor', [0 0.4470 0.7410], ...
            'MarkerSize', 8, 'LineWidth', 1.0, ...
            'DisplayName', 'Valid mass flow');
    end

    if ~isempty(invalidIdxPlot)
        plot(mDotList(invalidIdxPlot), (maxTcoldList(invalidIdxPlot) - 273.15), 'o', ...
            'MarkerEdgeColor', [0.6 0 0], ...
            'MarkerFaceColor', [0.85 0.33 0.10], ...
            'MarkerSize', 8, 'LineWidth', 1.0, ...
            'DisplayName', 'Invalid (T_{cold} > limit)');
    end

    xlabel('Water mass flow [kg/s]', ...
        'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Outlet water temperature [^\circC]', ...
        'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');

    formatAcademicAxes(gca);
    %lgd = legend('Location', 'best', 'Box', 'on', ...
        %'FontName', 'Times New Roman', 'FontSize', legendFontSize);
    %lgd.AutoUpdate = 'off';

    saveFigureAuto(fig, saveDir, 'water-outlet-temperature-vs-mass-flow');

    % If no valid reference exists, stop here
    if isempty(idxRef)
        return
    end

    xPlot      = x;
    TwaterPlot = T_CEA_ref(:,5) - 273.15;
    ThotPlot   = T_CEA_ref(:,2);
    TmidPlot   = T_CEA_ref(:,3);
    TcoldPlot  = T_CEA_ref(:,4);
    qPlot      = q_CEA_ref / 1e6;

    % 2) Hot wall temperature
    plotProfessionalVsX( ...
        xPlot, ThotPlot, ...
        'Axial coordinate, x [m]', 'Hot wall temperature [K]', ...
        'hot-wall-temperature', saveDir);

    % 3) TBC-metal interface temperature
    plotProfessionalVsX( ...
        xPlot, TmidPlot, ...
        'Axial coordinate, x [m]', 'TBC-metal interface temperature [K]', ...
        'tbc-metal-interface-temperature', saveDir);

    % 4) Coolant-side wall temperature
    plotProfessionalVsX( ...
        xPlot, TcoldPlot, ...
        'Axial coordinate, x [m]', 'Coolant-side wall temperature [K]', ...
        'coolant-side-wall-temperature', saveDir);

    % 5) Coolant temperature
    plotProfessionalVsX( ...
        xPlot, TwaterPlot, ...
        'Axial coordinate, x [m]', 'Coolant temperature ºC', ...
        'coolant-temperature-along-cooling-jacket', saveDir);

    % 6) Heat flux
    plotProfessionalVsX( ...
        xPlot, qPlot, ...
        'Axial coordinate, x [m]', 'Heat flux [MW/m^2]', ...
        'heat-flux-along-nozzle-axis', saveDir);

    % 7) All mass flows
    fig = figure('Color', 'w', 'Position', [100, 100, 850, 520], ...
             'Name', 'heat-flux-for-all-mass-flows');
    hold on;

    idxNonEmpty = find(~cellfun(@isempty, q_CEA_all));

    if ~isempty(idxNonEmpty)
        idxLowest  = idxNonEmpty(1);
        idxHighest = idxNonEmpty(end);
    else
        idxLowest  = [];
        idxHighest = [];
    end
    
    for j = 1:nCases
        if ~isempty(q_CEA_all{j})
            yData = q_CEA_all{j}/1e6; % Heat flux in MW/m^2
    
            % Determine label text for the lowest and highest mass flow
            if j == idxLowest
                labelTxt = sprintf('Mass flow = %.2f kg/s', mDotList(j));
                showInLegend = true;
            elseif j == idxHighest
                labelTxt = sprintf('Mass flow = %.2f kg/s', mDotList(j));
                showInLegend = true;
            else
                labelTxt = '';
                showInLegend = false;
            end
    
            % Plotting
            if validMdotMask(j)
                if showInLegend
                    plot(x, yData, '-', 'LineWidth', 1.8, ...
                        'DisplayName', labelTxt);
                else
                    plot(x, yData, '-', 'LineWidth', 1.4, 'HandleVisibility', 'off');
                end
            else
                if showInLegend
                    plot(x, yData, '--', 'LineWidth', 1.8, ...
                        'DisplayName', [labelTxt ' (invalid)']);
                else
                    plot(x, yData, '--', 'LineWidth', 1.4, 'HandleVisibility', 'off');
                end
            end
        end
end

xline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, 'HandleVisibility', 'off');

% Axis labels with LaTeX notation
xlabel('Axial coordinate, x [m]', ...
    'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Heat flux [MW/m^2]', ...
    'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');

formatAcademicAxes(gca);

% Legend configuration
lgd = legend('Location', 'east', 'Box', 'on', ...
    'FontName', 'Times New Roman', 'FontSize', legendFontSize);
lgd.AutoUpdate = 'off';

% Save the figure
saveFigureAuto(fig, saveDir, 'heat-flux-for-all-mass-flows');

% hola
%% ------------------------------------------------------------------------
% 1) Outlet water temperature vs mass flow (constant qDot = qDotMax)
% -------------------------------------------------------------------------
if design.found

    mFlow = 1:10;          % [kg/s]
    Tini_C = 18;           % [°C]
    ceWater = 4180;        % [J/kg/K]

    qDot = design.maxQDot; % [W/m^2]
    A = 61268.71789e-6;    % [m^2]

    Tfinal_C = Tini_C + (qDot * A) ./ (mFlow * ceWater);

    saveDir = 'C:\Users\celia\Desktop\figuras_cooling';
    legendFontSize = 18;

    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    fig = figure('Color', 'w', 'Position', [100, 100, 850, 520], ...
                 'Name', 'water-outlet-temperature-vs-mass-flow-design-MAX');
    hold on;

    % Plot for Outlet Water Temperature (constant qDot)
    plot(mFlow, Tfinal_C, '-', ...
        'Color', [0.4 0.4 0.4], 'LineWidth', 2, ...
        'HandleVisibility', 'off');  % Línea gris

    plot(mFlow, Tfinal_C, 'o', ...
        'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], ...  % Puntos azules
        'MarkerFaceColor', [0.8500, 0.3250, 0.0980], ...
        'MarkerSize', 8, 'LineWidth', 1.5, ...
        'DisplayName', 'Constant heat flux');  % Puntos azules

    % Plot for Real Outlet Temperature (from the CEA data)
    validIdxPlot   = find(validMdotMask);
    invalidIdxPlot = find(~validMdotMask);

    if ~isempty(validIdxPlot)
        plot(mDotList(validIdxPlot), TwaterOutList(validIdxPlot) - 273.15, '-', ...
            'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'HandleVisibility', 'off');  % Línea verde
        plot(mDotList(validIdxPlot), TwaterOutList(validIdxPlot) - 273.15, 'o', ...
            'MarkerEdgeColor', [0, 0.4470, 0.7410], ...  % Puntos verdes
            'MarkerFaceColor', [0, 0.4470, 0.7410], ...
            'MarkerSize', 8, 'LineWidth', 1.5, ...
            'DisplayName', 'Real mass flow');  % Puntos válidos
    end

    if ~isempty(invalidIdxPlot)
        plot(mDotList(invalidIdxPlot), (maxTcoldList(invalidIdxPlot) - 273.15), 'o', ...
            'MarkerEdgeColor', [0.8500 0.3250 0.0980], ...  % Puntos naranjas
            'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
            'MarkerSize', 8, 'LineWidth', 1.5, ...
            'DisplayName', 'Invalid (T_{cold} > limit)');  % Puntos inválidos
    end

    xlabel('Water mass flow [kg/s]', ...
        'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Outlet water temperature [^\circC]', ...
        'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');

    formatAcademicAxes(gca);

    % Añadimos la leyenda
    lgd = legend('Location', 'best', 'Box', 'on', ...
        'FontName', 'Times New Roman', 'FontSize', legendFontSize);
    lgd.AutoUpdate = 'off';

    % Guardar la figura
    saveFigureAuto(fig, saveDir, 'water-outlet-temperature-vs-mass-flow-design');
end
end
end

%% ========================================================================
% Solve one station iteratively because Bartz needs Thot
%% ========================================================================
function [q_i, hGas_i, Thot_i, Tmid_i, Tcold_i] = ...
    solveStationIterative(Taw_i, Tcool_i, ...
                          pc, cStar, Dt, rCurvature, epsilon_i, ...
                          mu_i, cp_i, Pr_i, ...
                          Tc_i, Te_i, ...
                          hWater_i, tWall, tTBC, kWall, kTBC)

Taw_i       = Taw_i(1);
Tcool_i     = Tcool_i(1);
pc          = pc(1);
cStar       = cStar(1);
Dt          = Dt(1);
rCurvature  = rCurvature(1);
epsilon_i   = epsilon_i(1);
mu_i        = mu_i(1);
cp_i        = cp_i(1);
Pr_i        = Pr_i(1);
Tc_i        = Tc_i(1);
Te_i        = Te_i(1);
hWater_i    = hWater_i(1);
tWall       = tWall(1);
tTBC        = tTBC(1);
kWall       = kWall(1);
kTBC        = kTBC(1);

Thot_i = 0.5*(Taw_i + Tcool_i);

for iter = 1:50
    hGas_i = bartzCorrelation(pc, cStar, Dt, rCurvature, epsilon_i, ...
                              mu_i, cp_i, Pr_i, Thot_i, Tc_i, Te_i, 0.6);

    hGas_i = hGas_i(1);

    H_i = 1 / (1/hGas_i + tWall/kWall + tTBC/kTBC + 1/hWater_i);
    q_i = H_i * (Taw_i - Tcool_i);

    Thot_new = Taw_i - q_i/hGas_i;

    if abs(Thot_new - Thot_i) < 1e-6
        Thot_i = Thot_new;
        break
    end

    Thot_i = 0.5*(Thot_i + Thot_new);
end

hGas_i = bartzCorrelation(pc, cStar, Dt, rCurvature, epsilon_i, ...
                          mu_i, cp_i, Pr_i, Thot_i, Tc_i, Te_i, 0.6);

hGas_i = hGas_i(1);

H_i = 1 / (1/hGas_i + tWall/kWall + tTBC/kTBC + 1/hWater_i);
q_i = H_i * (Taw_i - Tcool_i);

Thot_i = Taw_i - q_i/hGas_i;
Tmid_i = Thot_i - q_i*tTBC/kTBC;
Tcold_i = Tmid_i - q_i*tWall/kWall;

end


%% ========================================================================
% Very simple water-side HTC from fixed mass flow
%% ========================================================================
function [hWater, Re] = coolantHTCfromMassFlow(mDot, rhoWater, muWater, kWater, PrWater)

Dc = 1e-3;
Aflow = pi*(Dc^2)/4;
u = mDot/(rhoWater*Aflow);
Re = rhoWater*u*Dc/muWater;
Nu = 0.023*Re^0.8*PrWater^0.4;
hWater = Nu*kWater/Dc;

end


%% ========================================================================
% Simple CEA interpolation on mesh with reconstructed throat
%% ========================================================================
function [machCEA, TstatCEA, gammaCEA, cpCEA, muCEA, PrCEA] = ...
    buildCEAProfileOnMeshSimple(epsilon, n, pcNominal, muFallback, PrFallback)

iTh = n;

epsConvMesh = epsilon(1:iTh-1);
epsDivMesh  = epsilon(iTh+1:end);

epsConvCEA = sort(unique(epsConvMesh), 'descend').';
epsDivCEA  = sort(unique(epsDivMesh), 'ascend').';

[~, GasCEA] = getThermoProfileCEA_froz(80, 20, pcNominal/1e5, epsConvCEA, epsDivCEA);

nConv = numel(epsConvCEA);
nDiv  = numel(epsDivCEA);

epsConvData   = GasCEA.eps(1:nConv);
machConvData  = GasCEA.Mach(1:nConv);
TConvData     = GasCEA.T(1:nConv);
gammaConvData = GasCEA.gamma(1:nConv);
cpConvData    = GasCEA.cp(1:nConv);

epsDivData    = GasCEA.eps(nConv+1:nConv+nDiv);
machDivData   = GasCEA.Mach(nConv+1:nConv+nDiv);
TDivData      = GasCEA.T(nConv+1:nConv+nDiv);
gammaDivData  = GasCEA.gamma(nConv+1:nConv+nDiv);
cpDivData     = GasCEA.cp(nConv+1:nConv+nDiv);

machCEA  = zeros(size(epsilon));
TstatCEA = zeros(size(epsilon));
gammaCEA = zeros(size(epsilon));
cpCEA    = zeros(size(epsilon));

machCEA(1:iTh-1)  = interp1(epsConvData, machConvData, epsConvMesh, 'linear', 'extrap');
TstatCEA(1:iTh-1) = interp1(epsConvData, TConvData,    epsConvMesh, 'linear', 'extrap');
gammaCEA(1:iTh-1) = interp1(epsConvData, gammaConvData, epsConvMesh, 'linear', 'extrap');
cpCEA(1:iTh-1)    = interp1(epsConvData, cpConvData,   epsConvMesh, 'linear', 'extrap');

machCEA(iTh+1:end)  = interp1(epsDivData, machDivData, epsDivMesh, 'linear', 'extrap');
TstatCEA(iTh+1:end) = interp1(epsDivData, TDivData,    epsDivMesh, 'linear', 'extrap');
gammaCEA(iTh+1:end) = interp1(epsDivData, gammaDivData, epsDivMesh, 'linear', 'extrap');
cpCEA(iTh+1:end)    = interp1(epsDivData, cpDivData,   epsDivMesh, 'linear', 'extrap');

machCEA(iTh)  = 1.0;
TstatCEA(iTh) = 0.5*(TstatCEA(iTh-1) + TstatCEA(iTh+1));
gammaCEA(iTh) = 0.5*(gammaCEA(iTh-1) + gammaCEA(iTh+1));
cpCEA(iTh)    = 0.5*(cpCEA(iTh-1) + cpCEA(iTh+1));

if isfield(GasCEA,'viscosity')
    muConvData = GasCEA.viscosity(1:nConv);
    muDivData  = GasCEA.viscosity(nConv+1:nConv+nDiv);

    muCEA = zeros(size(epsilon));
    muCEA(1:iTh-1)   = interp1(epsConvData, muConvData, epsConvMesh, 'linear', 'extrap');
    muCEA(iTh+1:end) = interp1(epsDivData, muDivData,  epsDivMesh,  'linear', 'extrap');
    muCEA(iTh) = 0.5*(muCEA(iTh-1) + muCEA(iTh+1));
else
    muCEA = muFallback * ones(size(epsilon));
end

if isfield(GasCEA,'prandtl')
    PrConvData = GasCEA.prandtl(1:nConv);
    PrDivData  = GasCEA.prandtl(nConv+1:nConv+nDiv);

    PrCEA = zeros(size(epsilon));
    PrCEA(1:iTh-1)   = interp1(epsConvData, PrConvData, epsConvMesh, 'linear', 'extrap');
    PrCEA(iTh+1:end) = interp1(epsDivData, PrDivData,  epsDivMesh,  'linear', 'extrap');
    PrCEA(iTh) = 0.5*(PrCEA(iTh-1) + PrCEA(iTh+1));
else
    PrCEA = PrFallback * ones(size(epsilon));
end

end

%% ========================================================================
% Academic axis formatting
%% ========================================================================
function formatAcademicAxes(ax)

set(ax, 'FontName', 'Times New Roman', ...
        'FontSize', 12, ...
        'LineWidth', 0.8, ...
        'GridLineStyle', '--', ...
        'GridAlpha', 0.30, ...
        'Layer', 'top', ...
        'Box', 'on');

grid(ax, 'on');

end

%% ========================================================================
% Professional plot vs x
%% ========================================================================
function plotProfessionalVsX(xData, yData, xLabelTxt, yLabelTxt, fileName, saveDir)

fontName       = 'Times New Roman';
fontSize       = 12;
legendFontSize = 18;

techBlue   = [0, 0.4470, 0.7410];
techOrange = [0.8500, 0.3250, 0.0980];
lineColor  = [0.65, 0.65, 0.65];
markerSize = 7;
lineWidth  = 1.5;

fig = figure('Color', 'w', 'Position', [100, 100, 850, 520], 'Name', fileName);
hold on;

plot(xData, yData, '-', 'Color', lineColor, ...
     'LineWidth', lineWidth, 'HandleVisibility', 'off');

[~, idx_th] = min(abs(xData));
idx_left  = 1:(idx_th-1);
idx_right = (idx_th+1):length(xData);

if ~isempty(idx_left)
    h_left = plot(xData(idx_left), yData(idx_left), 'o', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'MarkerFaceColor', techBlue, ...
        'MarkerSize', markerSize, ...
        'LineWidth', 1, ...
        'DisplayName', 'Upstream stations');
else
    h_left = [];
end

h_th = plot(xData(idx_th), yData(idx_th), 'o', ...
    'MarkerEdgeColor', [0.2 0.2 0.2], ...
    'MarkerFaceColor', techOrange, ...
    'MarkerSize', markerSize + 2, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Throat');

if ~isempty(idx_right)
    h_right = plot(xData(idx_right), yData(idx_right), 'o', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'MarkerFaceColor', techBlue, ...
        'MarkerSize', markerSize, ...
        'LineWidth', 1, ...
        'DisplayName', 'Downstream stations');
else
    h_right = [];
end

xline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, 'HandleVisibility', 'off');

xlabel(xLabelTxt, 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
ylabel(yLabelTxt, 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');

xlim([min(xData), max(xData)]);

y_min = min(yData);
y_max = max(yData);
if y_max > y_min
    pad = 0.05*(y_max - y_min);
else
    pad = max(1, 0.05*abs(y_max));
end
ylim([y_min - pad, y_max + pad]);

formatAcademicAxes(gca);

if isempty(h_left)
    lgd = legend([h_th, h_right], 'Location', 'southeast', 'Box', 'on', ...
        'FontName', fontName, 'FontSize', legendFontSize);
elseif isempty(h_right)
    lgd = legend([h_left, h_th], 'Location', 'northeast', 'Box', 'on', ...
        'FontName', fontName, 'FontSize', legendFontSize);
else
    lgd = legend([h_left, h_th, h_right], 'Location', 'northeast', 'Box', 'on', ...
        'FontName', fontName, 'FontSize', legendFontSize);
end

lgd.AutoUpdate = 'off';

saveFigureAuto(fig, saveDir, fileName);

end
%% The png thingy
function saveFigureAuto(figHandle, saveDir, fileName)

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

pdfFile = fullfile(saveDir, [fileName '.pdf']);
pngFile = fullfile(saveDir, [fileName '.png']);

try
    % PDF vectorial
    exportgraphics(figHandle, pdfFile, 'ContentType', 'vector');
    
    % PNG opcional por si también lo quieres
    exportgraphics(figHandle, pngFile, 'Resolution', 300);

catch
    % Fallback para MATLAB antiguos
    print(figHandle, pdfFile, '-dpdf', '-bestfit');
    print(figHandle, pngFile, '-dpng', '-r300');
end

end