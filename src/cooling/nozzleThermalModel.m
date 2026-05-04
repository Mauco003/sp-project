function out = nozzleThermalModel(design, input, propellant, nozzle, performance, cooling, constants, options)
% nozzleThermalModel - Computes the thermal model of the nozzle, given a design and performance
%
% SYNTAX:
%  out = nozzleThermalModel(design, input, propellant, nozzle, performance, cooling, constants, options)
%
% INPUT:
%  design:          struct containing the design configuration and data
%  input:           struct containing the input parameters for the design
%  propellant:      struct containing the propellant properties
%  nozzle:          struct containing the nozzle geometry and properties
%  performance:     struct containing the performance data of the nozzle
%  cooling:         struct containing the cooling configuration and properties
%  constants:       struct containing the physical constants
%  options:         struct containing options for the thermal model computation:
%                       .n number of stations for discretization
%                       .areaLimits 2x1 vector with area limits for cooling
%                       .makePlot boolean to enable/disable plots
%                       .showSummary boolean to enable/disable summary print
%                       .savePlots boolean to enable/disable plot saving
%                       .mDotList vector of mass flow rates to evaluate
%                       .GasCEA optional struct with CEA data to bypass computation
%
% OUTPUT:
%  out: struct containing the results of the thermal model computation


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
        options.makePlot    (1,1) logical = false
        options.showSummary (1,1) logical = true
        options.savePlots   (1,1) logical = false
        options.mDotList    (1,:) double  = linspace(1,10,10)
        options.GasCEA                    = []
    end

    % INPUTS

    % from design

    tWall   = design.tWall;
    tTBC    = design.tTBC;
    kWall   = cooling(3).k;
    kTBC    = cooling(2).k;
    TwallMaxCoolantSide = cooling(4).Tboil;   % usually 150 C in K

    % properties of water

    ceWater  = cooling(4).cp;
    rhoWater = cooling(4).rho;
    muWater  = cooling(4).mu;
    kWater   = cooling(4).k;
    PrWater  = ceWater * muWater / kWater;

    % Inputs on Initial condition and mass flow

    waterT0       = cooling(4).Tini;
    waterPressure = cooling(4).pressure;

    mDotList = options.mDotList(:).';
    nCases   = numel(mDotList);

    % Chamber data

    R = constants.R / propellant.cea.molarMass;

    gamma = propellant.cea.gamma; % NOT constant along x
    pc    = input.pcNominal; % i need a fixed value for Pc, I choose nominal one - Steady state nozzle (time fixed)
    cStar = performance.cstar;
    muGas = propellant.cea.mu;
    cpGas = R * gamma/(gamma-1);
    PrGas = muGas*cpGas / propellant.cea.k;

    % Inputs from nozzle geometry

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

    % Heat exchange area of each station depends on the step on the x direction

    dA = 2*pi*(sqrt(diff(x).^2 + diff(r).^2)) .* 0.5 .* (r(1:end-1) + r(2:end));

    %% CEA MESH

    iTh = n;

    epsConvMesh = epsilon(1:iTh-1);
    epsDivMesh  = epsilon(iTh+1:end);

    epsConvCEA = sort(unique(epsConvMesh), 'descend').';
    epsDivCEA  = sort(unique(epsDivMesh), 'ascend').';


    if isempty(options.GasCEA)
        [~, GasCEA] = getThermoProfileCEAFroz(propellant.wtAP, propellant.wtHTPB, input.pcNominal/1e5, epsConvCEA, epsDivCEA);
    else
        GasCEA = options.GasCEA; % Bypass the expensive calls
    end

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

    % interpolate to save from singularities on the throat

    machCEA(1:iTh-1)  = interp1(epsConvData, machConvData, epsConvMesh, 'linear', 'extrap');
    TstatCEA(1:iTh-1) = interp1(epsConvData, TConvData, epsConvMesh, 'linear', 'extrap');
    gammaCEA(1:iTh-1) = interp1(epsConvData, gammaConvData, epsConvMesh, 'linear', 'extrap');
    cpCEA(1:iTh-1)    = interp1(epsConvData, cpConvData, epsConvMesh, 'linear', 'extrap');

    machCEA(iTh+1:end)  = interp1(epsDivData, machDivData, epsDivMesh, 'linear', 'extrap');
    TstatCEA(iTh+1:end) = interp1(epsDivData, TDivData, epsDivMesh, 'linear', 'extrap');
    gammaCEA(iTh+1:end) = interp1(epsDivData, gammaDivData, epsDivMesh, 'linear', 'extrap');
    cpCEA(iTh+1:end)    = interp1(epsDivData, cpDivData, epsDivMesh, 'linear', 'extrap');

    machCEA(iTh)  = 1.0;
    TstatCEA(iTh) = 0.5*(TstatCEA(iTh-1) + TstatCEA(iTh+1));
    gammaCEA(iTh) = 0.5*(gammaCEA(iTh-1) + gammaCEA(iTh+1));
    cpCEA(iTh)    = 0.5*(cpCEA(iTh-1) + cpCEA(iTh+1));

    if isfield(GasCEA,'viscosity')
        muConvData = GasCEA.viscosity(1:nConv);
        muDivData  = GasCEA.viscosity(nConv+1:nConv+nDiv);

        muCEA = zeros(size(epsilon));
        muCEA(1:iTh-1)   = interp1(epsConvData, muConvData, epsConvMesh, 'linear', 'extrap');
        muCEA(iTh+1:end) = interp1(epsDivData, muDivData, epsDivMesh, 'linear', 'extrap');
        muCEA(iTh) = 0.5*(muCEA(iTh-1) + muCEA(iTh+1));
    else
        muCEA = muGas * ones(size(epsilon));
    end

    if isfield(GasCEA,'prandtl')
        PrConvData = GasCEA.prandtl(1:nConv);
        PrDivData  = GasCEA.prandtl(nConv+1:nConv+nDiv);

        PrCEA = zeros(size(epsilon));
        PrCEA(1:iTh-1)   = interp1(epsConvData, PrConvData, epsConvMesh, 'linear', 'extrap');
        PrCEA(iTh+1:end) = interp1(epsDivData, PrDivData, epsDivMesh, 'linear', 'extrap');
        PrCEA(iTh) = 0.5*(PrCEA(iTh-1) + PrCEA(iTh+1));
    else
        PrCEA = PrGas * ones(size(epsilon));
    end

    % Try for different mass flows

    T_CEA_all        = cell(1,nCases);
    q_CEA_all        = cell(1,nCases);
    hGas_CEA_all     = cell(1,nCases);
    hWater_CEA_all   = cell(1,nCases);
    T0_CEA_all       = cell(1,nCases);
    maxQCEAList      = zeros(1,nCases);

    TwaterOutList    = NaN(1,nCases);
    maxTcoldList     = NaN(1,nCases);
    validMdotMask    = false(1,nCases);

    Dc = 1e-3;
    Aflow = pi*(Dc^2)/4;

    for j = 1:nCases

        mDot = mDotList(j);

        T_CEA      = zeros(N,5);   % [Taw, Thot, T_TBC_metal, Tcold, Twater]
        q_CEA      = zeros(N,1);
        hGas_CEA   = zeros(N,1);
        hWater_CEA = zeros(N,1);
        T0_CEA     = zeros(N,1);

        T_CEA(1,5) = waterT0;

        uWater = mDot/(rhoWater*Aflow);
        ReWater = rhoWater*uWater*Dc/muWater;
        hWater_i = dittusCorrelation(Dc,kWater,ReWater,PrWater);

        for i = 1:N

            T0_CEA(i)   = 2342.92;
            T_CEA(i,1)  = T0_CEA(i) * recoveryFactor(gammaCEA(i), machCEA(i), PrCEA(i));

            hWater_CEA(i) = hWater_i;

            Taw_i   = T_CEA(i,1);
            Tcool_i = T_CEA(i,5);

            Thot_i = 0.5*(Taw_i + Tcool_i);

            for iter = 1:50

                hGas_i = bartzCorrelation(pc, cStar, 2*rt, nozzle.rCurvature, epsilon(i), ...
                                        muCEA(i), cpCEA(i), PrCEA(i), ...
                                        Thot_i, T0_CEA(i), TstatCEA(i), 0.6);

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

            hGas_i = bartzCorrelation(pc, cStar, 2*rt, nozzle.rCurvature, epsilon(i), ...
                                    muCEA(i), cpCEA(i), PrCEA(i), ...
                                    Thot_i, T0_CEA(i), TstatCEA(i), 0.6);

            hGas_i = hGas_i(1);

            H_i = 1 / (1/hGas_i + tWall/kWall + tTBC/kTBC + 1/hWater_i);
            q_i = H_i * (Taw_i - Tcool_i);

            Thot_i = Taw_i - q_i/hGas_i;
            Tmid_i = Thot_i - q_i*tTBC/kTBC;
            Tcold_i = Tmid_i - q_i*tWall/kWall;

            q_CEA(i)    = q_i;
            hGas_CEA(i) = hGas_i;
            T_CEA(i,2)  = Thot_i;
            T_CEA(i,3)  = Tmid_i;
            T_CEA(i,4)  = Tcold_i;

            if i < N
                T_CEA(i+1,5) = T_CEA(i,5) + q_CEA(i)*dA(i)/(ceWater*mDot);
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

    %% PLOTS

    mDotRef = 2;

    idxRef = find(abs(mDotList - mDotRef) < 1e-12, 1);

    if ~isempty(idxRef)
        T_CEA_ref      = T_CEA_all{idxRef};
        q_CEA_ref      = q_CEA_all{idxRef};
        hGas_CEA_ref   = hGas_CEA_all{idxRef};
        hWater_CEA_ref = hWater_CEA_all{idxRef};
        T0_CEA_ref     = T0_CEA_all{idxRef};
    else
        T_CEA_ref      = [];
        q_CEA_ref      = [];
        hGas_CEA_ref   = [];
        hWater_CEA_ref = [];
        T0_CEA_ref     = [];
    end

    if ~isempty(q_CEA_ref)
        [maxQRef, idxMaxQRef] = max(q_CEA_ref);

        xTotal = max(x) - min(x);
        xOverL_maxQ = (x(idxMaxQRef) - min(x)) / xTotal;

        [~, idxThroat] = min(abs(x));
        xOverL_throat = (x(idxThroat) - min(x)) / xTotal;
    else
        maxQRef = NaN;
        idxMaxQRef = NaN;
        xOverL_maxQ = NaN;
        idxThroat = NaN;
        xOverL_throat = NaN;
    end

    %% OUTPUTS

    out = struct();

    % Add this to the output in order to run a single time the nozzleThemalModel in the main
    % and have the GasCEA without having to make the computation again everytime (if respected some conditions)
    out.GasCEA = GasCEA;

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
    out.hWaterCEA_all = hWater_CEA_all;
    out.T0_CEA_all = T0_CEA_all;

    out.machCEA   = machCEA;
    out.TstatCEA  = TstatCEA;
    out.gammaCEA  = gammaCEA;
    out.cpCEA     = cpCEA;
    out.muCEA     = muCEA;
    out.PrCEA     = PrCEA;

    out.maxQCEARef = maxQRef;
    out.idxMaxQCEARef = idxMaxQRef;
    out.xOverLMaxQCEARef = xOverL_maxQ;

    out.idxThroat = idxThroat;
    out.xOverLThroat = xOverL_throat;

    out.design = design;
    out.waterInletTemperature = waterT0;
    out.waterPressure = waterPressure;

    %% PRINTS

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
            fprintf('Reference max q               : %.4e W/m^2\n', maxQRef);
            fprintf('x/L at max q                  : %.4f\n', xOverL_maxQ);
            fprintf('x/L at throat                 : %.4f\n', xOverL_throat);
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

    %% PLOTS

    saveDir = fullfile(fileparts(mfilename("fullpath")), "..", "..", "figures");

    if options.makePlot

        legendFontSize = 18;

        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end

        % Outlet water temperature vs mass flow
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

        set(gca, 'FontName', 'Times New Roman', ...
                'FontSize', 12, ...
                'LineWidth', 0.8, ...
                'GridLineStyle', '--', ...
                'GridAlpha', 0.30, ...
                'Layer', 'top', ...
                'Box', 'on');
        grid on;

        if options.savePlots
            saveFigureAuto(fig, saveDir, 'water-outlet-temperature-vs-mass-flow');
        end

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

        plotList = { ...
            ThotPlot,   'Hot wall temperature [K]',           'hot-wall-temperature'; ...
            TmidPlot,   'TBC-metal interface temperature [K]', 'tbc-metal-interface-temperature'; ...
            TcoldPlot,  'Coolant-side wall temperature [K]',   'coolant-side-wall-temperature'; ...
            TwaterPlot, 'Coolant temperature ºC',              'coolant-temperature-along-cooling-jacket'; ...
            qPlot,      'Heat flux [MW/m^2]',                  'heat-flux-along-nozzle-axis'};

        for p = 1:size(plotList,1)

            yData = plotList{p,1};
            yLabelTxt = plotList{p,2};
            fileName = plotList{p,3};

            fig = figure('Color', 'w', 'Position', [100, 100, 850, 520], 'Name', fileName);
            hold on;

            plot(xPlot, yData, '-', 'Color', [0.65, 0.65, 0.65], ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');

            [~, idx_th] = min(abs(xPlot));
            idx_left  = 1:(idx_th-1);
            idx_right = (idx_th+1):length(xPlot);

            if ~isempty(idx_left)
                h_left = plot(xPlot(idx_left), yData(idx_left), 'o', ...
                    'MarkerEdgeColor', [0.2 0.2 0.2], ...
                    'MarkerFaceColor', [0, 0.4470, 0.7410], ...
                    'MarkerSize', 7, ...
                    'LineWidth', 1, ...
                    'DisplayName', 'Upstream stations');
            else
                h_left = [];
            end

            h_th = plot(xPlot(idx_th), yData(idx_th), 'o', ...
                'MarkerEdgeColor', [0.2 0.2 0.2], ...
                'MarkerFaceColor', [0.8500, 0.3250, 0.0980], ...
                'MarkerSize', 9, ...
                'LineWidth', 1.2, ...
                'DisplayName', 'Throat');

            if ~isempty(idx_right)
                h_right = plot(xPlot(idx_right), yData(idx_right), 'o', ...
                    'MarkerEdgeColor', [0.2 0.2 0.2], ...
                    'MarkerFaceColor', [0, 0.4470, 0.7410], ...
                    'MarkerSize', 7, ...
                    'LineWidth', 1, ...
                    'DisplayName', 'Downstream stations');
            else
                h_right = [];
            end

            xline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, 'HandleVisibility', 'off');

            xlabel('Axial coordinate, x [m]', ...
                'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel(yLabelTxt, ...
                'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');

            xlim([min(xPlot), max(xPlot)]);

            y_min = min(yData);
            y_max = max(yData);
            if y_max > y_min
                pad = 0.05*(y_max - y_min);
            else
                pad = max(1, 0.05*abs(y_max));
            end
            ylim([y_min - pad, y_max + pad]);

            set(gca, 'FontName', 'Times New Roman', ...
                    'FontSize', 12, ...
                    'LineWidth', 0.8, ...
                    'GridLineStyle', '--', ...
                    'GridAlpha', 0.30, ...
                    'Layer', 'top', ...
                    'Box', 'on');
            grid on;

            if isempty(h_left)
                lgd = legend([h_th, h_right], 'Location', 'southeast', 'Box', 'on', ...
                    'FontName', 'Times New Roman', 'FontSize', legendFontSize);
            elseif isempty(h_right)
                lgd = legend([h_left, h_th], 'Location', 'northeast', 'Box', 'on', ...
                    'FontName', 'Times New Roman', 'FontSize', legendFontSize);
            else
                lgd = legend([h_left, h_th, h_right], 'Location', 'northeast', 'Box', 'on', ...
                    'FontName', 'Times New Roman', 'FontSize', legendFontSize);
            end

            lgd.AutoUpdate = 'off';

            if options.savePlots
                saveFigureAuto(fig, saveDir, fileName);
            end
        end

        % All mass flows
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

        xlabel('Axial coordinate, x [m]', ...
            'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Heat flux [MW/m^2]', ...
            'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');

        set(gca, 'FontName', 'Times New Roman', ...
                'FontSize', 12, ...
                'LineWidth', 0.8, ...
                'GridLineStyle', '--', ...
                'GridAlpha', 0.30, ...
                'Layer', 'top', ...
                'Box', 'on');
        grid on;

        lgd = legend('Location', 'east', 'Box', 'on', ...
            'FontName', 'Times New Roman', 'FontSize', legendFontSize);
        lgd.AutoUpdate = 'off';

        if options.savePlots
            saveFigureAuto(fig, saveDir, 'heat-flux-for-all-mass-flows');
        end

        % hola
        %% ------------------------------------------------------------------------
        % 1) Outlet water temperature vs mass flow (constant qDot = qDotMax)
        % -------------------------------------------------------------------------
        if design.found

            mFlow = 1:10;          % [kg/s]

            qDot = design.maxQDot; % [W/m^2]

            % This area is hard coded because to take into accound the rounded throat it
            % has been obtained directly from a CAD software

            A_design = 61268.71789e-6;    % [m^2]

            waterTfKelvin = waterT0 + (qDot * A_design) ./ (mFlow * ceWater);
            waterTfCelsius = waterTfKelvin - 273.15;

            fig = figure('Color', 'w', 'Position', [100, 100, 850, 520], ...
                        'Name', 'water-outlet-temperature-vs-mass-flow-design-MAX');
            hold on;

            plot(mFlow, waterTfCelsius, '-', ...
                'Color', [0.4 0.4 0.4], 'LineWidth', 2, ...
                'HandleVisibility', 'off');

            plot(mFlow, waterTfCelsius, 'o', ...
                'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], ...
                'MarkerFaceColor', [0.8500, 0.3250, 0.0980], ...
                'MarkerSize', 8, 'LineWidth', 1.5, ...
                'DisplayName', 'Constant heat flux');

            validIdxPlot   = find(validMdotMask);
            invalidIdxPlot = find(~validMdotMask);

            if ~isempty(validIdxPlot)
                plot(mDotList(validIdxPlot), TwaterOutList(validIdxPlot) - 273.15, '-', ...
                    'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'HandleVisibility', 'off');
                plot(mDotList(validIdxPlot), TwaterOutList(validIdxPlot) - 273.15, 'o', ...
                    'MarkerEdgeColor', [0, 0.4470, 0.7410], ...
                    'MarkerFaceColor', [0, 0.4470, 0.7410], ...
                    'MarkerSize', 8, 'LineWidth', 1.5, ...
                    'DisplayName', 'Real mass flow');
            end

            if ~isempty(invalidIdxPlot)
                plot(mDotList(invalidIdxPlot), (maxTcoldList(invalidIdxPlot) - 273.15), 'o', ...
                    'MarkerEdgeColor', [0.8500 0.3250 0.0980], ...
                    'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
                    'MarkerSize', 8, 'LineWidth', 1.5, ...
                    'DisplayName', 'Invalid (T_{cold} > limit)');
            end

            xlabel('Water mass flow [kg/s]', ...
                'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Outlet water temperature [^\circC]', ...
                'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');

            set(gca, 'FontName', 'Times New Roman', ...
                    'FontSize', 12, ...
                    'LineWidth', 0.8, ...
                    'GridLineStyle', '--', ...
                    'GridAlpha', 0.30, ...
                    'Layer', 'top', ...
                    'Box', 'on');
            grid on;

            lgd = legend('Location', 'best', 'Box', 'on', ...
                'FontName', 'Times New Roman', 'FontSize', legendFontSize);
            lgd.AutoUpdate = 'off';

            if options.savePlots
                saveFigureAuto(fig, saveDir, 'water-outlet-temperature-vs-mass-flow-design');
            end
        end
    end

end
