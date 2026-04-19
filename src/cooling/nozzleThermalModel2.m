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
    options.compareCEA  (1,1) logical = true
end

%% ------------------------------------------------------------------------
% 0) Get fixed design values from design model
% -------------------------------------------------------------------------

tWall = design.tWall;
tTBC  = design.tTBC;

%% ------------------------------------------------------------------------
% 1) Fixed coolant data
% -------------------------------------------------------------------------
cpWater  = 4180;
rhoWater = 997;
muWater  = 0.89e-3;
kWater   = 0.60;
PrWater  = cpWater * muWater / kWater;

waterT0 = 18 + 273;
waterPressure = 10e5;
mDot = 3; 
%% ------------------------------------------------------------------------
% 2) Gas/chamber data
% -------------------------------------------------------------------------
R = constants.R / propellant.cea.molarMass;

gamma = propellant.cea.gamma;
Tc    = propellant.cea.ccTemperature;
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
if options.compareCEA
    [machCEA, TstatCEA, gammaCEA, cpCEA, muCEA, PrCEA] = ...
        buildCEAProfileOnMeshSimple(epsilon, n, input.pcNominal, muGas, PrGas);
else
    machCEA  = [];
    TstatCEA = [];
    gammaCEA = [];
    cpCEA    = [];
    muCEA    = [];
    PrCEA    = [];
end

%% ------------------------------------------------------------------------
% 5) Ideal/perfect-gas branch
% -------------------------------------------------------------------------
mach = zeros(N,1);
mach(n) = 1;
for i = 1:n-1
    mach(i) = machFromAreaRatio(epsilon(i), gamma, 'subsonic');
end
for i = n+1:N
    mach(i) = machFromAreaRatio(epsilon(i), gamma, 'supersonic');
end

T0 = 2342.92; % HARD CODED

% T(:,1)=Taw, T(:,2)=Thot, T(:,3)=after TBC, T(:,4)=cold metal, T(:,5)=Twater
T = zeros(N,5);
q = zeros(N,1);
hGas_ideal = zeros(N,1);
hWater_ideal = zeros(N,1);

T(1,5) = waterT0;

for i = 1:N
    Taw_i = T0 * recoveryFactor(gamma, mach(i), PrGas);
    T(i,1) = Taw_i;

    Tstat_i = T0 / (1 + 0.5*(gamma-1)*mach(i)^2);

    [hWater_i, ~] = coolantHTCfromMassFlow(mDot, rhoWater, muWater, kWater, PrWater);
    hWater_ideal(i) = hWater_i;

    [q_i, hGas_i, Thot_i, Tmid_i, Tcold_i] = ...
        solveStationIterative( ...
            Taw_i, T(i,5), ...
            pc, cStar, 2*rt, nozzle.rCurvature, epsilon(i), ...
            muGas, cpGas, PrGas, ...
            T0, Tstat_i, ...
            hWater_i, tWall, tTBC, design.kMetal, design.kTBC);

    q(i) = q_i;
    hGas_ideal(i) = hGas_i;
    T(i,2) = Thot_i;
    T(i,3) = Tmid_i;
    T(i,4) = Tcold_i;

    if i < N
        T(i+1,5) = T(i,5) + q(i)*dA(i)/(cpWater*mDot);
    end
end

%% ------------------------------------------------------------------------
% 6) CEA branch
% -------------------------------------------------------------------------
if options.compareCEA
    T_CEA  = zeros(N,5);
    q_CEA  = zeros(N,1);
    hGas_CEA = zeros(N,1);
    hWater_CEA = zeros(N,1);
    T0_CEA = zeros(N,1);

    T_CEA(1,5) = waterT0;

    for i = 1:N
        T0_CEA(i) = TstatCEA(i) * (1 + 0.5*(gammaCEA(i)-1)*machCEA(i)^2);
        T_CEA(i,1) = T0_CEA(i) * recoveryFactor(gammaCEA(i), machCEA(i), PrCEA(i));

        [hWater_i, ~] = coolantHTCfromMassFlow(mDot, rhoWater, muWater, kWater, PrWater);
        hWater_CEA(i) = hWater_i;

        [q_i, hGas_i, Thot_i, Tmid_i, Tcold_i] = ...
            solveStationIterative( ...
                T_CEA(i,1), T_CEA(i,5), ...
                input.pcNominal, performance.cstar, 2*rt, nozzle.rCurvature, epsilon(i), ...
                muCEA(i), cpCEA(i), PrCEA(i), ...
                T0_CEA(i), TstatCEA(i), ...
                hWater_i, tWall, tTBC, design.kMetal, design.kTBC);

        q_CEA(i) = q_i;
        hGas_CEA(i) = hGas_i;
        T_CEA(i,2) = Thot_i;
        T_CEA(i,3) = Tmid_i;
        T_CEA(i,4) = Tcold_i;

        if i < N
            T_CEA(i+1,5) = T_CEA(i,5) + q_CEA(i)*dA(i)/(cpWater*mDot);
        end
    end
else
    T_CEA  = [];
    q_CEA  = [];
    hGas_CEA = [];
    hWater_CEA = [];
    T0_CEA = [];
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
out.mDot  = mDot;

out.TIdeal = T;
out.qIdeal = q;
out.hGasIdeal = hGas_ideal;
out.hWaterIdeal = hWater_ideal;

out.TCEA = T_CEA;
out.qCEA = q_CEA;
out.hGasCEA = hGas_CEA;
out.hWaterCEA = hWater_CEA;
out.T0_CEA = T0_CEA;

out.machIdeal = mach;
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
    fprintf('=============================================\n');
    fprintf('      NOZZLE THERMAL MODEL - DISCRETIZED     \n');
    fprintf('=============================================\n');
    fprintf('tWall                 : %.6e m\n', tWall);
    fprintf('tTBC                  : %.6e m\n', tTBC);
    fprintf('mDot water            : %.4f kg/s\n', mDot);
    fprintf('Water inlet temp      : %.2f C\n', waterT0 - 273.15);
    fprintf('Max q ideal           : %.4e W/m^2\n', max(q));
    if options.compareCEA
        fprintf('Max q CEA             : %.4e W/m^2\n', max(q_CEA));
    end
    fprintf('=============================================\n\n');
end

%% ------------------------------------------------------------------------
% 9) Plots
% -------------------------------------------------------------------------
if options.makePlot

    if options.compareCEA
        xPlot      = x;
        TwaterPlot = T_CEA(:,5) - 273.15;
        TawPlot    = T_CEA(:,1);
        ThotPlot   = T_CEA(:,2);
        TmidPlot   = T_CEA(:,3);
        TcoldPlot  = T_CEA(:,4);
        qPlot      = q_CEA / 1e6;     % MW/m^2
        hGasPlot   = hGas_CEA;
        machPlot   = machCEA;
        T0Plot     = T0_CEA;
        TstatPlot  = TstatCEA;
        gammaPlot  = gammaCEA;
        cpPlot     = cpCEA;
    else
        xPlot      = x;
        TwaterPlot = T(:,5) - 273.15;
        TawPlot    = T(:,1);
        ThotPlot   = T(:,2);
        TmidPlot   = T(:,3);
        TcoldPlot  = T(:,4);
        qPlot      = q / 1e6;         % MW/m^2
        hGasPlot   = hGas_ideal;
        machPlot   = mach;
        T0Plot     = [];
        TstatPlot  = [];
        gammaPlot  = [];
        cpPlot     = [];
    end

    % 1) Coolant temperature
    plotProfessionalVsX( ...
        xPlot, TwaterPlot, ...
        'Axial coordinate, x [m]', 'Coolant Temperature [^\circC]', ...
        'Coolant temperature along cooling jacket', ...
        'Coolant');

    % 2) Hot wall temperature
    plotProfessionalVsX( ...
        xPlot, ThotPlot, ...
        'Axial coordinate, x [m]', 'Hot Wall Temperature [K]', ...
        'Hot-side wall temperature', ...
        'Hot Wall');

    % 3) Temperatures through layers
    f = figure('Color', 'w', 'Position', [100, 100, 900, 540], ...
               'Name', 'Temperature through layers');
    hold on;

    p1 = plot(xPlot, TawPlot,   '-',  'LineWidth', 1.8, 'DisplayName', 'T_{aw}');
    p2 = plot(xPlot, ThotPlot,  '-',  'LineWidth', 1.8, 'DisplayName', 'T_{hot}');
    p3 = plot(xPlot, TmidPlot,  '-',  'LineWidth', 1.8, 'DisplayName', 'T_{TBC-metal}');
    p4 = plot(xPlot, TcoldPlot, '-',  'LineWidth', 1.8, 'DisplayName', 'T_{cold metal}');
    p5 = plot(xPlot, TwaterPlot + 273.15, '-', 'LineWidth', 1.8, 'DisplayName', 'T_{water}');
    xline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, 'HandleVisibility', 'off');

    xlabel('Axial coordinate, x [m]', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Temperature [K]',         'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
    title('Temperature distribution across gas, wall and coolant', ...
          'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');

    formatAcademicAxes(gca);
    legend([p1 p2 p3 p4 p5], 'Location', 'eastoutside', 'Box', 'on', 'FontName', 'Times New Roman');

    % 4) Heat flux
    plotProfessionalVsX( ...
        xPlot, qPlot, ...
        'Axial coordinate, x [m]', 'q'''' [MW/m^2]', ...
        'Heat flux along nozzle axis', ...
        'Heat Flux');

    % 5) Gas-side HTC
    plotProfessionalVsX( ...
        xPlot, hGasPlot, ...
        'Axial coordinate, x [m]', 'h_g [W/m^2/K]', ...
        'Gas-side heat transfer coefficient along nozzle', ...
        'h_g');

    % 6) Mach number
    plotProfessionalVsX( ...
        xPlot, machPlot, ...
        'Axial coordinate, x [m]', 'Mach [-]', ...
        'Mach number along nozzle', ...
        'Mach');

    % 7) CEA stagnation and static temperatures
    if ~isempty(T0Plot)
        f = figure('Color', 'w', 'Position', [100, 100, 850, 520], ...
                   'Name', 'CEA temperatures');
        hold on;

        p1 = plot(xPlot, T0Plot,    '-', 'LineWidth', 1.8, 'DisplayName', 'T_0');
        p2 = plot(xPlot, TstatPlot, '-', 'LineWidth', 1.8, 'DisplayName', 'T_{stat}');
        xline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, 'HandleVisibility', 'off');

        xlabel('Axial coordinate, x [m]', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Temperature [K]',         'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
        title('CEA stagnation and static temperatures', ...
              'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');

        formatAcademicAxes(gca);
        legend([p1 p2], 'Location', 'east', 'Box', 'on', 'FontName', 'Times New Roman');
    end

    % 8) Gamma
    if ~isempty(gammaPlot)
        plotProfessionalVsX( ...
            xPlot, gammaPlot, ...
            'Axial coordinate, x [m]', '\gamma [-]', ...
            'CEA heat capacity ratio', ...
            '\gamma');
    end

    % 9) cp
    if ~isempty(cpPlot)
        plotProfessionalVsX( ...
            xPlot, cpPlot, ...
            'Axial coordinate, x [m]', 'c_p [J/kg/K]', ...
            'CEA specific heat at constant pressure', ...
            'c_p');
    end

    % 10) Final qdot plot requested
    plotProfessionalVsX( ...
        xPlot, qPlot, ...
        'Axial coordinate, x [m]', 'q'''' [MW/m^2]', ...
        'Local heat flux distribution q''''(x)', ...
        'q''''');

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

% Force everything to be scalar
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

    hGas_i = hGas_i(1);   % in case Bartz returns a 1x1 or vector

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
Tmid_i = Thot_i - q_i*tWall/kWall;
Tcold_i = Tmid_i - q_i*tTBC/kTBC;

end


%% ========================================================================
% Very simple water-side HTC from fixed mass flow
%% ========================================================================
function [hWater, Re] = coolantHTCfromMassFlow(mDot, rhoWater, muWater, kWater, PrWater)

Dc = 1e-3; % fixed simple characteristic diameter
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

[~, GasCEA] = getThermoProfileCEA(80, 20, pcNominal/1e5, epsConvCEA, epsDivCEA);

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
function plotProfessionalVsX(xData, yData, xLabelTxt, yLabelTxt, titleTxt, displayNameTxt)

fontName   = 'Times New Roman';
fontSize   = 12;
techBlue   = [0, 0.4470, 0.7410];
lineColor  = [0.65, 0.65, 0.65];
markerSize = 7;
lineWidth  = 1.5;

f = figure('Color', 'w', 'Position', [100, 100, 850, 520], 'Name', titleTxt);
hold on;

% Grey path
plot(xData, yData, '-', 'Color', lineColor, ...
     'LineWidth', lineWidth, 'HandleVisibility', 'off');

% Throat index
[~, idx_th] = min(abs(xData));

idx_left  = 1:(idx_th-1);
idx_right = (idx_th+1):length(xData);

% Left branch
if ~isempty(idx_left)
    h_left = plot(xData(idx_left), yData(idx_left), 's', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'MarkerFaceColor', techBlue, ...
        'MarkerSize', markerSize, ...
        'LineWidth', 1, ...
        'DisplayName', 'Upstream Stations');
else
    h_left = [];
end

% Throat
techOrange = [0.8500, 0.3250, 0.0980];
h_th = plot(xData(idx_th), yData(idx_th), 'o', ...
    'MarkerEdgeColor', [0.2 0.2 0.2], ...
    'MarkerFaceColor', techOrange, ...
    'MarkerSize', markerSize + 2, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Throat');

% Right branch
if ~isempty(idx_right)
    h_right = plot(xData(idx_right), yData(idx_right), 'd', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'MarkerFaceColor', techBlue, ...
        'MarkerSize', markerSize, ...
        'LineWidth', 1, ...
        'DisplayName', 'Downstream Stations');
else
    h_right = [];
end

xline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, 'HandleVisibility', 'off');

xlabel(xLabelTxt, 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
ylabel(yLabelTxt, 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
title(titleTxt, 'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');

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
    legend([h_th, h_right], 'Location', 'east', 'Box', 'on', 'FontName', fontName);
elseif isempty(h_right)
    legend([h_left, h_th], 'Location', 'east', 'Box', 'on', 'FontName', fontName);
else
    legend([h_left, h_th, h_right], 'Location', 'east', 'Box', 'on', 'FontName', fontName);
end

end
