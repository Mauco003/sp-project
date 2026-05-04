function [] = plotPerformance(performance, performanceCEA, viewSavedPlots, options)
% plotPerformance - Plots the performance parameters of the rocket motor over time, comparing simplified calculations with CEA results.
%
% SYNTAX:
%  plotPerformance(performance, performanceCEA, viewSavedPlots, options)
%
% INPUT:
%  performance - A structure containing the simplified performance parameters over time.
%  performanceCEA - A structure containing the CEA performance parameters over time.
%  viewSavedPlots - A string ('on' or 'off') to control the visibility of saved plots.
%  options.savePlots - A logical flag to control whether to save the generated plots as PDF files.

    arguments
        performance, performanceCEA, 
        viewSavedPlots = 'off'
        options.savePlots (1,1) logical = true
    end
    
    % Extract time vector and ensure all data is in column format for consistent plotting
    t = performance.t(:); 
    N = length(t);
    cleanData = @(data) (isscalar(data) * ones(N, 1) + ~isscalar(data) * data(:));
    
    % Simplified parameters
    thrustSimp = cleanData(performance.thrust);
    ctSimp     = cleanData(performance.ct);
    cstarSimp  = cleanData(performance.cstar);
    mDotSimp   = cleanData(performance.mDot);
    IspSimp    = cleanData(performance.Isp);
    pcSimp     = cleanData(performance.pc); 
    peSimp     = cleanData(performance.pe);
    
    % CEA parameters
    thrustCEA = cleanData(performanceCEA.thrust);
    ctCEA     = cleanData(performanceCEA.ct);
    cstarCEA  = cleanData(performanceCEA.cstar);
    IspCEA    = cleanData(performanceCEA.Isp);
    peCEA     = cleanData(performanceCEA.pe);
    pcCEA     = cleanData(performanceCEA.pc);

    % Style settings
    fontName   = 'Times New Roman';
    fontSize   = 24;
    colorSim   = [0.40, 0.40, 0.40];       % Dark gray
    colorCEA   = [0, 0.4470, 0.7410];      % CEA technic blue
    colorAtm   = [0.34, 0.34, 0.34];       % Light gray for p_atm
    colorPe    = [0.8500, 0.3250, 0.0980]; % CEA technic orange for Pe
    lineWidth  = 1.5;
    
    % General plot (subplot 2x2)
    f_combined = figure('Color', 'w', 'Position', [50, 50, 1200, 900], 'Name', 'Performance Comparison'); 
    
    subplot(2, 2, 1); hold on; box on; grid on;
    plot(t, IspSimp, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, IspCEA, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Specific Impulse, I_{sp} [s]', fontName, fontSize, t, true);
    ylim([0, max([IspSimp; IspCEA]) * 1.05]);
    
    subplot(2, 2, 2); hold on; box on; grid on;
    plot(t, peSimp ./ 1e5, '--', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, peCEA ./ 1e5, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    yline(1.01325, ':', 'Color', colorAtm, 'LineWidth', lineWidth + 0.5, 'DisplayName', 'Pamb');
    formatAxis('Time [s]', 'Exit Pressure, Pe [bar]', fontName, fontSize, t, true);
    ylim([0, max([peSimp; peCEA; 1.01325*1e5] ./ 1e5) * 1.05]);
    
    subplot(2, 2, 3); hold on; box on; grid on;
    plot(t, pcSimp ./ 1e5, '--', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, pcCEA ./ 1e5, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Chamber Pressure, Pc [bar]', fontName, fontSize, t, true);
    ylim([0, max([pcSimp; pcCEA] ./ 1e5) * 1.05]);
    
    subplot(2, 2, 4); hold on; box on; grid on;
    plot(t, thrustSimp ./ 1000, '--', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, thrustCEA ./ 1000, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Thrust, T [kN]', fontName, fontSize, t, true);
    ylim([0, max([thrustSimp; thrustCEA] ./ 1000) * 1.05]);
    
    if options.savePlots
        exportgraphics(f_combined, 'Performance_Comparison_Combined.pdf', 'ContentType', 'vector');
    end

    % Save single plots

    % COMBINED PLOT Pc and Pe (YYAXIS) with 70:1 Alignment
    f_dual = figure('Color', 'w', 'Position', [150, 150, 900, 600], 'Visible', viewSavedPlots);
    hold on; 
    
    % Define the alignment target
    pc_target = 70;
    pe_target = 1;
    scaling_factor = pc_target / pe_target; % Results in 70
    
    % Calculate the common upper limit (based on the maximum value between the two)
    % Add a 10% margin for aesthetics
    max_pe_needed = max([peSimp; peCEA; 1.01325*1e5] ./ 1e5);
    max_pc_needed = max([pcSimp; pcCEA] ./ 1e5);
    
    % Choose the limit that satisfies both while maintaining the 70:1 ratio
    top_pe = max(max_pe_needed, max_pc_needed / scaling_factor) * 1.1;
    top_pc = top_pe * scaling_factor;

    % LEFT AXIS: Chamber Pressure
    yyaxis left
    plot(t, pcSimp ./ 1e5, '--', 'Color', colorSim, 'LineWidth', lineWidth, 'HandleVisibility', 'off');
    plot(t, pcCEA ./ 1e5, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'Pc CEA');
    ylabel('Chamber Pressure, Pc [bar]', 'FontName', fontName, 'FontSize', fontSize + 4, 'FontWeight', 'bold');
    set(gca, 'YColor', [0 0 0]);
    ylim([0, top_pc]); % Scaled limit

    % RIGHT AXIS: Exit Pressure
    yyaxis right
    plot(t, peSimp ./ 1e5, '--', 'Color', colorPe, 'LineWidth', lineWidth, 'DisplayName', 'Pe Sim');
    plot(t, peCEA ./ 1e5, '-', 'Color', colorPe, 'LineWidth', lineWidth, 'DisplayName', 'Pe CEA');
    yline(1.01325, ':', 'Color', colorAtm, 'LineWidth', lineWidth, 'DisplayName', 'Pamb');
    ylabel('Exit Pressure, Pe [bar]', 'FontName', fontName, 'FontSize', fontSize + 4, 'FontWeight', 'bold');
    set(gca, 'YColor', colorPe);
    ylim([0, top_pe]); % Scaled limit

    % Specific manual formatting
    ax = gca;
    ax.Box = 'on'; 
    grid on;
    set(ax, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 0.8, 'GridLineStyle', '--', 'GridAlpha', 0.3);
    xlabel('Time [s]', 'FontName', fontName, 'FontSize', fontSize + 4, 'FontWeight', 'bold');
    xlim([min(t), max(t)]);
    
    legend('Location', 'best', 'Box', 'on', 'FontSize', fontSize-6, 'FontName', fontName);
    if options.savePlots
        exportgraphics(f_dual, 'single_pc_pe_dual.pdf', 'ContentType', 'vector');
    end

    % Remaining Single Plots
    % Thrust (T)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, thrustSimp ./ 1000, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, thrustCEA ./ 1000, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Thrust, T [kN]', fontName, fontSize, t, true);
    ylim([0, max([thrustSimp; thrustCEA] ./ 1000) * 1.05]);
    if options.savePlots
        exportgraphics(f, 'single_thrust.pdf', 'ContentType', 'vector');
    end

    % Specific Impulse (Isp)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, IspSimp, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, IspCEA, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Specific Impulse, I_{sp} [s]', fontName, fontSize, t, true);
    ylim([0, max([IspSimp; IspCEA]) * 1.05]);
    if options.savePlots
        exportgraphics(f, 'single_isp.pdf', 'ContentType', 'vector');
    end

    % Thrust Coefficient (ct)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, ctSimp, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, ctCEA, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Thrust Coefficient [-]', fontName, fontSize, t, true);
    if options.savePlots
        exportgraphics(f, 'single_ct.pdf', 'ContentType', 'vector');
    end

    % Characteristic Velocity (c*)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, cstarSimp, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cstarCEA, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'c^* [m/s]', fontName, fontSize, t, true);
    if options.savePlots
        exportgraphics(f, 'single_cstar.pdf', 'ContentType', 'vector');
    end

    % Mass Flow Rate (mDot)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', viewSavedPlots); 
    hold on; box on; grid on;
    plot(t, mDotSimp, '-', 'Color', colorCEA, 'LineWidth', lineWidth);
    formatAxis('Time [s]', 'Mass Flow Rate [kg/s]', fontName, fontSize, t, false);
    ylim([0, max(mDotSimp) * 1.05]);
    if options.savePlots
        exportgraphics(f, 'single_mass_flow.pdf', 'ContentType', 'vector');
    end
end

function formatAxis(xLabelStr, yLabelStr, fontName, fontSize, t, showLegend)
    ax = gca;
    set(ax, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 0.8, ...
            'GridLineStyle', '--', 'GridAlpha', 0.3, 'Layer', 'top');
    labelSize = fontSize + 4;
    xlabel(xLabelStr, 'FontName', fontName, 'FontSize', labelSize, 'FontWeight', 'bold');
    ylabel(yLabelStr, 'FontName', fontName, 'FontSize', labelSize, 'FontWeight', 'bold');
    xlim([min(t), max(t)]);
    if showLegend
        legend('Location', 'best', 'Box', 'on', 'FontName', fontName, 'FontSize', fontSize-4); 
    end
    
end