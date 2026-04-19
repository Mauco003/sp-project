function [] = plotPerformance(performance, performanceCEA, viewSavedPlots)
    arguments
        performance, performanceCEA, 
        viewSavedPlots = 'off'
    end
    
    % --- 0. ALLINEAMENTO DIMENSIONALE DEI DATI ---
    t = performance.t(:); 
    N = length(t);
    cleanData = @(data) (isscalar(data) * ones(N, 1) + ~isscalar(data) * data(:));

    % Parametri Semplificati
    sim_thrust = cleanData(performance.thrust);
    sim_ct     = cleanData(performance.ct);
    sim_cstar  = cleanData(performance.cstar);
    sim_mDot   = cleanData(performance.mDot);
    sim_Isp    = cleanData(performance.Isp);
    sim_pc     = cleanData(performance.pc); 
    
    % Parametri CEA
    cea_thrust = cleanData(performanceCEA.thrust);
    cea_ct     = cleanData(performanceCEA.ct);
    cea_cstar  = cleanData(performanceCEA.cstar);
    cea_Isp    = cleanData(performanceCEA.Isp);
    cea_pe     = cleanData(performanceCEA.pe);
    cea_pc     = cleanData(performanceCEA.pc);

    % --- Style Settings ---
    fontName   = 'Times New Roman';
    fontSize   = 12;
    colorSim   = [0.40, 0.40, 0.40];       % Grigio scuro
    colorCEA   = [0, 0.4470, 0.7410];      % Blu tecnico CEA
    colorAtm   = [0.34, 0.34, 0.34];       % Grigio chiaro per p_atm
    lineWidth  = 1.5;
    
    % --- 1. GRAFICO GENERALE (SUBPLOTS 2x2) ---
    f_combined = figure('Color', 'w', 'Position', [50, 50, 950, 750], 'Name', 'Performance Comparison'); 
    
    % Top-Left: Specific Impulse (Isp) -> Y-LIM MODIFICATO
    subplot(2, 2, 1); hold on; box on; grid on;
    plot(t, sim_Isp, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_Isp, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Specific Impulse, I_{sp} [s]', fontName, fontSize, t, true);
    ylim([0, max([sim_Isp; cea_Isp]) * 1.05]);
    
    % Top-Right: Exit Pressure (pe) -> Y-LIM MODIFICATO
    subplot(2, 2, 2); hold on; box on; grid on;
    plot(t, cea_pe ./ 1e5, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    yline(1.01325, '--', 'Color', colorAtm, 'LineWidth', lineWidth + 0.5, 'DisplayName', 'P_{amb}');
    formatAxis('Time [s]', 'Exit Pressure, P_e [bar]', fontName, fontSize, t, true);
    legend('Location', 'southeast', 'Box', 'on', 'FontName', fontName);
    ylim([0, max([cea_pe ./ 1e5; 1.01325]) * 1.05]);

    % Bottom-Left: Chamber Pressure (pc) -> Y-LIM MODIFICATO
    subplot(2, 2, 3); hold on; box on; grid on;
    plot(t, sim_pc ./ 1e5, '--', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_pc ./ 1e5, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Chamber Pressure, P_c [bar]', fontName, fontSize, t, true);
    ylim([0, max([sim_pc; cea_pc] ./ 1e5) * 1.05]);

    % Bottom-Right: Thrust (T) -> Y-LIM MODIFICATO
    subplot(2, 2, 4); hold on; box on; grid on;
    plot(t, sim_thrust ./ 1000, '--', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_thrust ./ 1000, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Thrust, T [kN]', fontName, fontSize, t, true);
    ylim([0, max([sim_thrust; cea_thrust] ./ 1000) * 1.05]);

    exportgraphics(f_combined, 'Performance_Comparison_Combined.pdf', 'ContentType', 'vector');

    % --- 2. SALVATAGGIO GRAFICI SINGOLI (UNPACKED) ---

    % Thrust (T) -> Y-LIM MODIFICATO
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, sim_thrust ./ 1000, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_thrust ./ 1000, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Thrust, T [kN]', fontName, fontSize, t, true);
    ylim([0, max([sim_thrust; cea_thrust] ./ 1000) * 1.05]);
    exportgraphics(f, 'single_thrust.pdf', 'ContentType', 'vector');

    % Specific Impulse (Isp) -> Y-LIM MODIFICATO
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, sim_Isp, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_Isp, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Specific Impulse, I_{sp} [s]', fontName, fontSize, t, true);
    ylim([0, max([sim_Isp; cea_Isp]) * 1.05]);
    exportgraphics(f, 'single_isp.pdf', 'ContentType', 'vector');

    % Exit Pressure (Pe) -> Y-LIM MODIFICATO
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, cea_pe ./ 1e5, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    yline(1.01325, '-', 'Color', colorAtm, 'LineWidth', lineWidth + 0.5, 'DisplayName', 'P_{amb}');
    formatAxis('Time [s]', 'Exit Pressure, P_e [bar]', fontName, fontSize, t, true);
    legend('Location', 'west', 'Box', 'on', 'FontName', fontName);
    ylim([0, max([cea_pe ./ 1e5; 1.01325]) * 1.05]);
    exportgraphics(f, 'single_pe.pdf', 'ContentType', 'vector');

    % Chamber Pressure (Pc) -> Y-LIM MODIFICATO
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, sim_pc ./ 1e5, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_pc ./ 1e5, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Chamber Pressure, P_c [bar]', fontName, fontSize, t, true);
    ylim([0, max([sim_pc; cea_pc] ./ 1e5) * 1.05]);
    exportgraphics(f, 'single_pc.pdf', 'ContentType', 'vector');

    % --- GRAFICI AGGIUNTIVI (LIMITI STANDARD) ---

    % Thrust Coefficient (ct)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, sim_ct, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_ct, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Thrust Coefficient [-]', fontName, fontSize, t, true);
    exportgraphics(f, 'single_ct.pdf', 'ContentType', 'vector');

    % Characteristic Velocity (c*)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, sim_cstar, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_cstar, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'c^* [m/s]', fontName, fontSize, t, true);
    exportgraphics(f, 'single_cstar.pdf', 'ContentType', 'vector');

    % Mass Flow Rate (mDot)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, sim_mDot, '-', 'Color', colorCEA, 'LineWidth', lineWidth);
    formatAxis('Time [s]', 'Mass Flow Rate [kg/s]', fontName, fontSize, t, false);
    ylim([0, max(sim_mDot) * 1.05]);
    exportgraphics(f, 'single_mass_flow.pdf', 'ContentType', 'vector');
end

% --- Helper: Formattazione Assi ---
function formatAxis(xLabelStr, yLabelStr, fontName, fontSize, t, showLegend)
    xlabel(xLabelStr, 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
    ylabel(yLabelStr, 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
    xlim([min(t), max(t)]);
    ax = gca;
    set(ax, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 0.8, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'Layer', 'top');
    if showLegend, legend('Location', 'best', 'Box', 'on', 'FontName', fontName); end
end