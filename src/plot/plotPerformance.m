function [] = plotPerformance(performance, performanceCEA)
    
    % 0. ALLINEAMENTO DIMENSIONALE DEI DATI
    
    t = performance.t(:); 
    N = length(t);
    
    function v = cleanData(data)
        if isscalar(data)
            v = data * ones(N, 1);
        else
            v = data(:); 
        end
    end

    % Parametri Semplificati
    sim_thrust = cleanData(performance.thrust);
    sim_ct     = cleanData(performance.ct);
    sim_cstar  = cleanData(performance.cstar);
    sim_mDot   = cleanData(performance.mDot);
    sim_Isp    = cleanData(performance.Isp);

    % Parametri CEA
    cea_thrust = cleanData(performanceCEA.thrust);
    cea_ct     = cleanData(performanceCEA.ct);
    cea_cstar  = cleanData(performanceCEA.cstar);
    cea_Isp    = cleanData(performanceCEA.Isp);
    cea_pe     = cleanData(performanceCEA.pe);

    % --- Style Settings ---
    fontName   = 'Times New Roman';
    fontSize   = 12;
    colorSim   = [0.40, 0.40, 0.40];       % Grigio scuro
    colorCEA   = [0, 0.4470, 0.7410];      % Blu tecnico CEA
    colorAtm   = [0.34, 0.34, 0.34];       % Grigio chiaro per p_atm
    lineWidth  = 1.5;

    
    % 1. GRAFICO GENERALE (SUBPLOTS 2x2)
    
    f_combined = figure('Color', 'w', 'Position', [50, 50, 950, 750], 'Name', 'Performance Comparison'); 
    
    % --- Alto a sinistra: Specific Impulse (Isp) ---
    subplot(2, 2, 1); hold on; box on; grid on;
    plot(t, sim_Isp, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_Isp, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Specific Impulse, I_{sp} [s]', fontName, fontSize, t, true);
    
     % --- Alto a destra: Exit Pressure (pe) ---
    subplot(2, 2, 2); hold on; box on; grid on;
    plot(t, cea_pe ./ 1e5, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    yline(1.01325, '--', 'Color', colorAtm, 'LineWidth', lineWidth + 0.5, 'DisplayName', 'P_{amb}');
    formatAxis('Time [s]', 'Exit Pressure, P_e [bar]', fontName, fontSize, t, true);
    legend('Location', 'southeast', 'Box', 'on', 'FontName', fontName); % Forza la legenda in basso

    % --- Basso a sinistra: Characteristic Velocity (c*) ---
    subplot(2, 2, 3); hold on; box on; grid on;
    plot(t, sim_cstar, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_cstar, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'c^* [m/s]', fontName, fontSize, t, true);

    % --- Basso a destra: Thrust Coefficient (ct) ---
    subplot(2, 2, 4); hold on; box on; grid on;
    plot(t, sim_ct, '-', 'Color', colorSim, 'LineWidth', lineWidth, 'DisplayName', 'Simplified');
    plot(t, cea_ct, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    formatAxis('Time [s]', 'Thrust Coefficient [-]', fontName, fontSize, t, true);

    % Salvataggio combinato
    fileNameCombined = 'Performance_Comparison_Combined.pdf';
    try
        exportgraphics(f_combined, fileNameCombined, 'ContentType', 'vector');
    catch
        print(f_combined, fileNameCombined, '-dpdf', '-r300', '-bestfit');
    end

    
    % 2. SALVATAGGIO GRAFICI SINGOLI
  
    % Spinta (T)
    saveSinglePlot(t, sim_thrust ./ 1000, cea_thrust ./ 1000, ...
        'Simplified', 'CEA', 'Time [s]', 'Thrust, T [kN]', 'Perf_Single_Thrust.pdf', ...
        colorSim, colorCEA, lineWidth, fontName, fontSize);

    % Coefficienti e Indici
    saveSinglePlot(t, sim_ct, cea_ct, 'Simplified', 'CEA', 'Time [s]', 'Thrust Coefficient [-]', 'Perf_Single_ThrustCoeff.pdf', colorSim, colorCEA, lineWidth, fontName, fontSize);
    saveSinglePlot(t, sim_cstar, cea_cstar, 'Simplified', 'CEA', 'Time [s]', 'c^* [m/s]', 'Perf_Single_Cstar.pdf', colorSim, colorCEA, lineWidth, fontName, fontSize);
    saveSinglePlot(t, sim_Isp, cea_Isp, 'Simplified', 'CEA', 'Time [s]', 'Specific Impulse, I_{sp} [s]', 'Perf_Single_Isp.pdf', colorSim, colorCEA, lineWidth, fontName, fontSize);
    
    % Pressione di Uscita (Grafico speciale solo con CEA e p_atm)
    saveSinglePlot_Pe(t, cea_pe ./ 1e5, 'Time [s]', 'Exit Pressure, P_e [bar]', 'Perf_Single_Pe.pdf', colorCEA, colorAtm, lineWidth, fontName, fontSize);

    % Massa (mDot)
    saveSinglePlot_SingleLine_NoLegend(t, sim_mDot, 'Time [s]', 'Mass Flow Rate [kg/s]', 'Perf_Single_MassFlow.pdf', colorCEA, lineWidth, fontName, fontSize);
end

% --- Helper: Formattazione Assi ---
function formatAxis(xLabelStr, yLabelStr, fontName, fontSize, t, showLegend)
    xlabel(xLabelStr, 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
    ylabel(yLabelStr, 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
    xlim([min(t), max(t)]); % Limite tagliato esattamente alla fine del burn
    ax = gca;
    set(ax, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 0.8, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'Layer', 'top');
    if showLegend
        legend('Location', 'best', 'Box', 'on', 'FontName', fontName); 
    end
end

% --- Helper: Grafici Doppi ---
function saveSinglePlot(t, data1, data2, name1, name2, xLabelStr, yLabelStr, fileName, color1, color2, lineWidth, fontName, fontSize)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', 'off'); 
    hold on; box on; grid on;
    plot(t, data1, '-', 'Color', color1, 'LineWidth', lineWidth, 'DisplayName', name1);
    plot(t, data2, '-', 'Color', color2, 'LineWidth', lineWidth, 'DisplayName', name2);
    formatAxis(xLabelStr, yLabelStr, fontName, fontSize, t, true);
    exportgraphics(f, fileName, 'ContentType', 'vector');
    close(f);
end

% --- Helper: Grafico Singolo Pressione Uscita (CEA + p_atm) ---
function saveSinglePlot_Pe(t, dataCEA, xLabelStr, yLabelStr, fileName, colorCEA, colorAtm, lineWidth, fontName, fontSize)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', 'off'); 
    hold on; box on; grid on;
    
    plot(t, dataCEA, '-', 'Color', colorCEA, 'LineWidth', lineWidth, 'DisplayName', 'CEA');
    
    % Aumento lo spessore per far risaltare il tratteggio
    yline(1.01325, '-', 'Color', colorAtm, 'LineWidth', lineWidth + 0.5, 'DisplayName', 'P_{amb}');
    
    formatAxis(xLabelStr, yLabelStr, fontName, fontSize, t, true);
    
    % Forza la legenda in basso a destra
legend('Location', 'west', 'Box', 'on', 'FontName', fontName);    
    exportgraphics(f, fileName, 'ContentType', 'vector');
    close(f);
end

% --- Helper: Grafico Singolo mDot ---
function saveSinglePlot_SingleLine_NoLegend(t, data1, xLabelStr, yLabelStr, fileName, color1, lineWidth, fontName, fontSize)
    f = figure('Color', 'w', 'Position', [150, 150, 800, 400], 'Visible', 'off'); 
    hold on; box on; grid on;
    plot(t, data1, '-', 'Color', color1, 'LineWidth', lineWidth);
    formatAxis(xLabelStr, yLabelStr, fontName, fontSize, t, false);
    exportgraphics(f, fileName, 'ContentType', 'vector');
    close(f);
end