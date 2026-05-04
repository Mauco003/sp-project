function CJdesign = thermalModelDesign(input, nozzle, performance, cooling, options)
% thermalModelDesign - Design of the thermal protection system for the nozzle
%
% SYNTAX:
%  CJdesign = thermalModelDesign(input, nozzle, performance, cooling, options)
%
% INPUT:
%  input        - struct with design input parameters
%  nozzle       - struct with nozzle design parameters
%  performance  - struct with performance parameters from CEA
%  cooling      - struct with cooling configuration and material properties
%  options      - struct with optional parameters:
%                   .showPlots (logical) - whether to show plots (default: false)
%                   .showSummary (logical) - whether to show summary (default: true)
%                   .savePlots (logical) - whether to save plots (default: false)
%
% OUTPUT:
%  CJdesign     - struct with the results of the thermal design:
%                   .found (logical) - whether a valid design was found
%                   .THot (double) - design hot wall temperature [K]
%                   .tTBC (double) - required TBC thickness [m]
%                   .tWall (double) - required wall thickness [m]
%                   .criticalPhase (int) - index of the critical station
%                   .criticalName (string) - name of the critical station
%                   .qDotCritical (double) - heat flux at critical station [W/m^2]
%                   .THotCritical (double) - hot wall temperature at critical station [K]
%                   .THotMetalCritical (double) - metal temperature at critical station [K]
%                   .TColdMetalCritical (double) - coolant side metal temperature at critical station [K]
%                   .maxQDot (double) - maximum heat flux across all stations [W/m^2]
%                   .maxTHot (double) - maximum hot wall temperature across all stations [K]
%                   .minTHot (double) - minimum hot wall temperature across all stations [K]
%                   .maxTHotMetal (double) - maximum metal temperature across all stations [K]
%                   .minTHotMetal (double) - minimum metal temperature across all stations [K]



    arguments
        input struct
        nozzle struct
        performance struct
        cooling struct
        options.showPlots (1,1) logical = false
        options.showSummary (1,1) logical = true
        options.savePlots (1,1) logical = false
    end

    %% Get inputs

    % parameters obtained from nozzle design. In cooling there is the fixed
    % information on the material chosen for TBC, nozzle wall and also the
    % constants for water and hot flow obtained in CEA

    cfg.pc = input.pcNominal; % check this suspicious thing         
    cfg.T0 = 2342.92;              

    cfg.At = nozzle.At;
    cfg.Dt = 2*sqrt(cfg.At/pi);

    cfg.rCurvature = nozzle.rCurvature; % needed for Bartz correlation
    cfg.cStar = performance.cstar(1);

    % wall thickness calculation

    inconelYS = cooling(3).YS;
    safetyFactor = cooling.SF;
    cfg.tWall = 65e5*sqrt(2*cfg.At/pi)/inconelYS*safetyFactor;              

    cfg.kTBC   = cooling(2).k;               
    cfg.kMetal = cooling(3).k;              

    cfg.Thm_max = cooling(3).Tmax; % boiling of inconel  
    cfg.Tcw_max = cooling(4).Tboil; % boiling point of water 

    % maximum and minimum operation temperature

    cfg.THot_min = 1000; % realistic value           
    cfg.THot_max = cooling(2).Tmax;           
    cfg.dTHot    = 10;             

    % maximum and minimum thickness of TBC

    cfg.tTBC_min = 100e-6; % realistic value       
    cfg.tTBC_max = cooling(2).tmax;         

    % because I need to check the validity of the design on several stations
    % (because i don't know yet where the maximum heat flux will be and maybe
    % my design doesnt work in other stations, I assume it will be the throat
    % but just in case I check in more areas) 

    stationEps   = [2, 1.5, 1, 1.5, 2]; % stations to check ( 1 is the throat and will probably be the most critical)
    stationNames = {'Inlet e=2','Conv e=1.5','Throat','Div e=1.5','Exit e=2'};

    % Build the 5 stations from CEA

    [~, GasCEA_2]  = getThermoProfileCEAFroz(80, 20, cfg.pc/1e5, 2.0, 2.0);
    [~, GasCEA_15] = getThermoProfileCEAFroz(80, 20, cfg.pc/1e5, 1.5, 1.5);

    gammaVals = [GasCEA_2.gamma(1),     GasCEA_15.gamma(1),     GasCEA_2.gamma(2),     GasCEA_15.gamma(3),     GasCEA_2.gamma(3)];
    TVals     = [GasCEA_2.T(1),         GasCEA_15.T(1),         GasCEA_2.T(2),         GasCEA_15.T(3),         GasCEA_2.T(3)];
    muVals    = [GasCEA_2.viscosity(1), GasCEA_15.viscosity(1), GasCEA_2.viscosity(2), GasCEA_15.viscosity(3), GasCEA_2.viscosity(3)];
    cpVals    = [GasCEA_2.cp(1),        GasCEA_15.cp(1),        GasCEA_2.cp(2),        GasCEA_15.cp(3),        GasCEA_2.cp(3)];
    prVals    = [GasCEA_2.prandtl(1),   GasCEA_15.prandtl(1),   GasCEA_2.prandtl(2),   GasCEA_15.prandtl(3),   GasCEA_2.prandtl(3)];
    machVals  = [GasCEA_2.Mach(1),      GasCEA_15.Mach(1),      GasCEA_2.Mach(2),      GasCEA_15.Mach(3),      GasCEA_2.Mach(3)];

    n = length(stationEps);

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
        stations(i).epsilon   = stationEps(i);
        stations(i).gamma     = gammaVals(i);
        stations(i).Tstat     = TVals(i);
        stations(i).viscosity = muVals(i);
        stations(i).cp        = cpVals(i);
        stations(i).prandtl   = prVals(i);
        stations(i).Mach      = machVals(i);
    end


    % Calculation of thickness of TBC and operating temperature (hot wall)

    CJdesign.found = false;

    CJdesign.lastFail.THot = NaN;
    CJdesign.lastFail.requiredThickness = NaN;
    CJdesign.lastFail.station = NaN;
    CJdesign.lastFail.reason = '';

    nStations = length(stations);

    for THot_design = cfg.THot_min : cfg.dTHot : cfg.THot_max

        for i = 1:nStations

            station = stations(i);

            gamma = station.gamma;
            Tstat = station.Tstat;
            muGas = station.viscosity;
            cpGas = station.cp;
            PrGas = station.prandtl;
            mach  = station.Mach;
            epsi  = station.epsilon;
            
            % Get adiabatic wall temperature from hot flow temp from CEA

            r = recoveryFactor(gamma, mach, PrGas);
            Taw = cfg.T0*r;
            
            % Bartz Correlation using the fixed Thotwall (design parameter)

            hGas = bartzCorrelation( ...
                cfg.pc, cfg.cStar, cfg.Dt, ...
                cfg.rCurvature, epsi, ...
                muGas, cpGas, PrGas, ...
                THot_design, cfg.T0, Tstat, 0.6);
            
            % Convective heat transfer from hot flow to wall
            qDot = hGas*(Taw - THot_design);
            
            % q dot is kept constant (steady assumption) get the thickness
            % needed from the TBC for that q dot

            Rcond_req = (THot_design - cfg.Tcw_max)/qDot;
            tReq = cfg.kTBC*(Rcond_req - cfg.tWall/cfg.kMetal);
            
            % With given ddot and thickness of both TBC and wall, the
            % temperature on the hot side of the cooling duct is obtained. This
            % is our restricting parameter as a hot point on this wall would
            % make the water boil

            THotMetal = cfg.Tcw_max + qDot*cfg.tWall/cfg.kMetal;
            
            % TBC THICKNESS CONSTRAIN. It has to be greater than the minimum
            % value

            if tReq < cfg.tTBC_min

                CJdesign.lastFail.THot = THot_design;
                CJdesign.lastFail.requiredThickness = tReq;
                CJdesign.lastFail.station = i;
                CJdesign.lastFail.reason = ...
                    'Required thickness is smaller than tTBC_min';

                break

            end
            
            % TBC THICKNESS CONSTRAIN (physical limits given by material)

            if tReq > cfg.tTBC_max

                CJdesign.lastFail.THot = THot_design;
                CJdesign.lastFail.requiredThickness = tReq;
                CJdesign.lastFail.station = i;
                CJdesign.lastFail.reason = ...
                    'Required thickness is larger than tTBC_max';

                break

            end
            
            % Boiling water constraint (Hot wall of coolant side has a maximum
            % temperature)

            if THotMetal > cfg.Thm_max

                CJdesign.lastFail.THot = THot_design;
                CJdesign.lastFail.requiredThickness = tReq;
                CJdesign.lastFail.station = i;
                CJdesign.lastFail.reason = ...
                    'THot_target matched but THotMetal exceeds Thm_max';

                break

            end

            % check next station. It is assumed that the critical point will be
            % in the throat, but if we check multiple stations we can validate
            % this hipotesis. If the TBC works in the throat it will work
            % anywhere else

            % We repeat the HT with a fixed thickness. If it fails, then we increase the
            % ThotWall (design parameter)

            if i < nStations

                nextStation = stations(i+1);

                gamma = nextStation.gamma;
                Tstat = nextStation.Tstat;
                muGas = nextStation.viscosity;
                cpGas = nextStation.cp;
                PrGas = nextStation.prandtl;
                mach  = nextStation.Mach;
                epsi  = nextStation.epsilon;

                r = recoveryFactor(gamma, mach, PrGas);
                Taw = cfg.T0*r;

                Rcond = tReq/cfg.kTBC + cfg.tWall/cfg.kMetal;

                THot_guess = min(max(THot_design, cfg.Tcw_max + 1e-6), Taw - 1e-6);

                for iter = 1:50

                    hGas = bartzCorrelation( ...
                        cfg.pc, cfg.cStar, cfg.Dt, ...
                        cfg.rCurvature, epsi, ...
                        muGas, cpGas, PrGas, ...
                        THot_guess, cfg.T0, Tstat, 0.6);

                    qDot_next = (Taw - cfg.Tcw_max)/(1/hGas + Rcond);

                    nextTHot = cfg.Tcw_max + qDot_next*Rcond;
                    nextTHotMetal = cfg.Tcw_max + qDot_next*cfg.tWall/cfg.kMetal;

                    if abs(nextTHot - THot_guess) < 1e-8
                        break
                    end

                    THot_guess = 0.5*THot_guess + 0.5*nextTHot;

                end

                nextWorks = ...
                    (nextTHot <= THot_design + 1e-6) && ...
                    (nextTHotMetal <= cfg.Thm_max + 1e-6);

            else

                nextWorks = true;

            end

            if nextWorks

                CJdesign.found = true;

                CJdesign.THot = THot_design;
                CJdesign.tTBC = tReq;
                CJdesign.tWall = cfg.tWall;

                CJdesign.criticalPhase = i;
                CJdesign.criticalName = station.name;

                CJdesign.qDotCritical = qDot;
                CJdesign.THotCritical = THot_design;
                CJdesign.THotMetalCritical = THotMetal;
                CJdesign.TColdMetalCritical = cfg.Tcw_max;

                CJdesign.maxQDot = qDot;

                CJdesign.maxTHot = THot_design;
                CJdesign.minTHot = THot_design;

                CJdesign.maxTHotMetal = THotMetal;
                CJdesign.minTHotMetal = THotMetal;

                CJdesign.maxTColdMetal = cfg.Tcw_max;
                CJdesign.minTColdMetal = cfg.Tcw_max;

                break

            end

        end

        if CJdesign.found
            break
        end

    end


    % When the design has been found, the trade off between hot wall, thickness
    % and maximum coolant wall hot side temperature has been found

    % Assuming that the maximum heat flux is mantained constant along the
    % jacket I can estimamte the final temperature outside of the cooling
    % jacket with Q = mass*ce*A(Tini_water - Tfin_water) by trying with
    % different mass flows 

    if options.showPlots
    if CJdesign.found

        mFlow = 1:10;          % [kg/s]
        Tini_C = 18;           % [°C]
        ceWater = 4180;        % [J/kg/K]

        qDot = CJdesign.maxQDot; % [W/m^2]
        A = 61268.71789e-6;    % [m^2]

        Tfinal_C = Tini_C + (qDot * A) ./ (mFlow * ceWater);

        saveDir = fullfile(fileparts(mfilename("fullpath")), "..", "..", "figures");

        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end

        fig = figure('Color', 'w', 'Position', [100, 100, 850, 520], ...
                    'Name', 'water-outlet-temperature-vs-mass-flow-design-MAX');
        hold on;

        plot(mFlow, Tfinal_C, '-', ...
            'Color', [0.65 0.65 0.65], 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');

        plot(mFlow, Tfinal_C, 'o', ...
            'MarkerEdgeColor', [0.2 0.2 0.2], ...
            'MarkerFaceColor', [0 0.4470 0.7410], ...
            'MarkerSize', 8, 'LineWidth', 1.0, ...
            'DisplayName', 'Constant q'''' = q''''_{max}');

        xlabel('Water mass flow [kg/s]', ...
            'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Outlet water temperature [^\circC]', ...
            'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');

        formatAcademicAxes(gca);

        if options.savePlots
            saveFigureAuto(fig, saveDir, 'water-outlet-temperature-vs-mass-flow-design');
        end
    end
    end

    if options.showSummary
        printSummary(CJdesign);
    end
end

% Print and save results

function printSummary(best)

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

end



%% Helper functions
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