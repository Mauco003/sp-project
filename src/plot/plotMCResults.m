function plotMCResults(mcData, config)
% plotMCResults - Generates all visualizations for the Monte Carlo analysis.
% This function strictly handles UI/plotting and performs no physics computations.
%
% SYNTAX:
%  plotMCResults(mcData, config)
%
% INPUT:
%  mcData - A structure containing all Monte Carlo results and inputs
%  config - A structure containing configuration settings for plotting and analysis modes


    %% Unpack Configuration for cleaner code
    fontName   = config.plots.fontName;
    fontSize   = config.plots.fontSize;
    techBlue   = config.plots.techBlue;
    techOrange = config.plots.techOrange;
    techGreen  = config.plots.techGreen;
    techPurple = config.plots.techPurple;
    lineWidth  = config.plots.lineWidth;
    
    totalRunsMC = mcData.totalRunsMC;
    analysisMode = config.analysisMode;
    
    %% Detailed Convergence Monitoring
    validBT = ~isnan(mcData.burnTime);
    validMEOP = ~isnan(mcData.MEOP);
    
    burnTimeMeanProg = cumsum(mcData.burnTime, 'omitnan') ./ cumsum(validBT);
    MEOPMeanProg = cumsum(mcData.MEOP, 'omitnan') ./ cumsum(validMEOP);
    
    burnTimeStdProg = zeros(totalRunsMC, 1);
    MEOPStdProg = zeros(totalRunsMC, 1);
    for k = 2:totalRunsMC
        burnTimeStdProg(k) = std(mcData.burnTime(1:k), 'omitnan');
        MEOPStdProg(k) = std(mcData.MEOP(1:k), 'omitnan');
    end
    
    tolPercent = 5; 
    limitLabel = sprintf('%g%% LIMIT', tolPercent);
    tolValueBTMean = mcData.nominals.burningTime * (tolPercent / 100);
    tolValueBTStd  = mcData.burnTimeStd * (tolPercent / 100);
    tolValueMEOPMean = (mcData.nominals.pcNominal * 1e-5) * (tolPercent / 100);
    tolValueMEOPStd  = (mcData.MEOPStd * 1e-5) * (tolPercent / 100);
    
    figure('Color', 'w', 'Name', 'Detailed Convergence Monitoring', 'WindowState', 'maximized');
    
    % Subplot 1: BT MEAN
    subplot(2, 2, 1); hold on; grid on; box on;
    plot(1:totalRunsMC, burnTimeMeanProg, '-', 'Color', techBlue, 'LineWidth', lineWidth);
    yline(mcData.nominals.burningTime, '-k', 'DESIGN', 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'bottom', 'LineWidth', lineWidth);
    yline(mcData.nominals.burningTime + tolValueBTMean, '-', limitLabel, 'Color', techOrange, 'LabelHorizontalAlignment', 'center', 'LineWidth', lineWidth);
    yline(mcData.nominals.burningTime - tolValueBTMean, '-', 'Color', techOrange, 'LineWidth', lineWidth);
    xlabel('Iterations [-]'); ylabel('BT MEAN [s]'); xlim([1, totalRunsMC]);
    yMin = min([mcData.nominals.burningTime - 1.1*tolValueBTMean, min(burnTimeMeanProg, [], 'omitnan')]);
    yMax = max([mcData.nominals.burningTime + 1.1*tolValueBTMean, max(burnTimeMeanProg, [], 'omitnan')]);
    ylim([yMin, yMax]);
    set(gca, 'FontName', fontName, 'FontSize', fontSize, 'GridAlpha', 0.3);
    
    % Subplot 2: BT STD
    subplot(2, 2, 2); hold on; grid on; box on;
    plot(1:totalRunsMC, burnTimeStdProg, '-', 'Color', techBlue, 'LineWidth', lineWidth);
    yline(mcData.burnTimeStd, '-k', 'FINAL STD', 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'bottom', 'LineWidth', lineWidth);
    yline(mcData.burnTimeStd + tolValueBTStd, '-', limitLabel, 'Color', techOrange, 'LabelHorizontalAlignment', 'center', 'LineWidth', lineWidth);
    yline(mcData.burnTimeStd - tolValueBTStd, '-', 'Color', techOrange, 'LineWidth', lineWidth);
    xlabel('Iterations [-]'); ylabel('BT STD [s]'); xlim([1, totalRunsMC]);
    yMin = max(0, min([mcData.burnTimeStd - 1.1*tolValueBTStd, min(burnTimeStdProg, [], 'omitnan')]));
    yMax = max([mcData.burnTimeStd + 1.1*tolValueBTStd, max(burnTimeStdProg, [], 'omitnan')]);
    ylim([yMin, yMax]);
    set(gca, 'FontName', fontName, 'FontSize', fontSize, 'GridAlpha', 0.3);
    
    % Subplot 3: MEOP MEAN
    subplot(2, 2, 3); hold on; grid on; box on;
    plot(1:totalRunsMC, MEOPMeanProg*1e-5, '-', 'Color', techBlue, 'LineWidth', lineWidth);
    yline(mcData.nominals.pcNominal*1e-5, '-k', 'DESIGN', 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'bottom', 'LineWidth', lineWidth);
    yline(mcData.nominals.pcNominal*1e-5 + tolValueMEOPMean, '-', limitLabel, 'Color', techOrange, 'LabelHorizontalAlignment', 'center', 'LineWidth', lineWidth);
    yline(mcData.nominals.pcNominal*1e-5 - tolValueMEOPMean, '-', 'Color', techOrange, 'LineWidth', lineWidth);
    xlabel('Iterations [-]'); ylabel('MEOP MEAN [bar]'); xlim([1, totalRunsMC]);
    yMin = min([(mcData.nominals.pcNominal*1e-5) - 1.1*tolValueMEOPMean, min(MEOPMeanProg*1e-5, [], 'omitnan')]);
    yMax = max([(mcData.nominals.pcNominal*1e-5) + 1.1*tolValueMEOPMean, max(MEOPMeanProg*1e-5, [], 'omitnan')]);
    ylim([yMin, yMax]);
    set(gca, 'FontName', fontName, 'FontSize', fontSize, 'GridAlpha', 0.3);
    
    % Subplot 4: MEOP STD
    subplot(2, 2, 4); hold on; grid on; box on;
    plot(1:totalRunsMC, MEOPStdProg*1e-5, '-', 'Color', techBlue, 'LineWidth', lineWidth);
    yline(mcData.MEOPStd*1e-5, '-k', 'FINAL STD', 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'bottom', 'LineWidth', lineWidth);
    yline(mcData.MEOPStd*1e-5 + tolValueMEOPStd, '-', limitLabel, 'Color', techOrange, 'LabelHorizontalAlignment', 'center', 'LineWidth', lineWidth);
    yline(mcData.MEOPStd*1e-5 - tolValueMEOPStd, '-', 'Color', techOrange, 'LineWidth', lineWidth);
    xlabel('Iterations [-]'); ylabel('MEOP STD [bar]'); xlim([1, totalRunsMC]);
    yMin = max(0, min([(mcData.MEOPStd*1e-5) - 1.1*tolValueMEOPStd, min(MEOPStdProg*1e-5, [], 'omitnan')]));
    yMax = max([(mcData.MEOPStd*1e-5) + 1.1*tolValueMEOPStd, max(MEOPStdProg*1e-5, [], 'omitnan')]);
    ylim([yMin, yMax]);
    set(gca, 'FontName', fontName, 'FontSize', fontSize, 'GridAlpha', 0.3);

    %% Primary Histograms
    figure('Color', 'w', 'Name', 'Primary Output Distributions');
    
    subplot(1, 2, 1); hold on; grid on; box on;
    histogram(mcData.MEOP*1e-5, 'Normalization', 'pdf', 'FaceColor', techGreen, 'EdgeColor', 'w');
    plotNormalFit(mcData.MEOP*1e-5); 
    xline(mcData.MEOPMean*1e-5, '--', 'MEAN', 'Color', techOrange, 'LineWidth', lineWidth);
    xlabel('MEOP [bar]', 'FontWeight', 'bold'); ylabel('PDF [-]', 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', fontSize, 'GridAlpha', 0.3);
    
    subplot(1, 2, 2); hold on; grid on; box on;
    histogram(mcData.burnTime, 'Normalization', 'pdf', 'FaceColor', techPurple, 'EdgeColor', 'w');
    plotNormalFit(mcData.burnTime); 
    xline(mcData.burnTimeMean, '--', 'MEAN', 'Color', techOrange, 'LineWidth', lineWidth);
    xlabel('Burning Time [s]', 'FontWeight', 'bold'); ylabel('PDF [-]', 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', fontSize, 'GridAlpha', 0.3);

    %% 2D Scatter Landscapes
    if analysisMode == 0 || analysisMode == 1
        figure('Color', 'w', 'Name', '2D Design Space Scatter', 'WindowState', 'maximized');
        
        subplot(1, 2, 1); hold on; box on; grid on;
        scatter(mcData.inputs.a, mcData.inputs.n, 40, mcData.MEOP * 1e-5, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        cb1 = colorbar; cb1.Label.String = 'MEOP [bar]'; colormap(gca, jet); 
        xline(mcData.nominals.a, '--', 'Nominal ''a'''); yline(mcData.nominals.n, '--', 'Nominal ''n''');
        xlabel('Burn Rate Coefficient (a)'); ylabel('Burn Rate Exponent (n)'); title('Pressure Landscape');
        set(gca, 'FontName', fontName, 'FontSize', 14, 'GridAlpha', 0.3);
    
        subplot(1, 2, 2); hold on; box on; grid on;
        scatter(mcData.inputs.a, mcData.inputs.n, 40, mcData.burnTime, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        cb2 = colorbar; cb2.Label.String = 'Burn Time [s]'; colormap(gca, flipud(jet)); 
        xline(mcData.nominals.a, '--', 'Nominal ''a'''); yline(mcData.nominals.n, '--', 'Nominal ''n''');
        xlabel('Burn Rate Coefficient (a)'); ylabel('Burn Rate Exponent (n)'); title('Burn Time Landscape');
        set(gca, 'FontName', fontName, 'FontSize', 14, 'GridAlpha', 0.3);
    end

    %% Comprehensive Plots & Sensitivities
    if config.useComprehensive
        % Extended Distributions
        figure('Color','w','Name','Extended Output Distributions');
        
        subplot(2,2,1); hold on;
        plotRobustHistogram(mcData.TMax, techBlue, 'T_{max}', fontName);
        
        subplot(2,2,2); hold on;
        plotRobustHistogram(mcData.IspAvg, techOrange, 'Isp_{avg}', fontName);
        
        subplot(2,2,3); hold on;
        plotRobustHistogram(mcData.IspTotMax, techGreen, 'Total Impulse', fontName);
        
        subplot(2,2,4); hold on;
        plotRobustHistogram(mcData.thrustReliability, techPurple, 'Reliability (>90% T) [%]', fontName);
        
        sgtitle('Extended Monte Carlo Outputs', 'FontName', fontName, 'FontWeight', 'bold');
        
        % Relative Boxplot
        devMEOP = (mcData.MEOP - mcData.MEOPMean) ./ mcData.MEOPMean * 100;
        devTMax = (mcData.TMax - mean(mcData.TMax)) ./ mean(mcData.TMax) * 100;
        devIspAvg = (mcData.IspAvg - mean(mcData.IspAvg)) ./ mean(mcData.IspAvg) * 100;
        devBurnTime = (mcData.burnTime - mcData.burnTimeMean) ./ mcData.burnTimeMean * 100;

        figure('Color', 'w', 'Name', 'Normalized Monte Carlo Spread');
        boxplot([devMEOP(:), devTMax(:), devIspAvg(:), devBurnTime(:)], ...
            'Labels', {'MEOP', 'T_{max}', 'Isp_{avg}', 'Burn Time'}, 'Colors', techBlue, 'Symbol', 'rx');
        ylabel('Deviation from Mean [%]', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
        title('Relative Performance Uncertainty', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
        grid on; set(gca, 'FontName', fontName, 'FontSize', fontSize, 'GridLineStyle', '--', 'GridAlpha', 0.3);

        % Input Matrix Generation for Sensitivity
        inputMatrix = [mcData.inputs.a, mcData.inputs.n, mcData.inputs.OF, mcData.inputs.pAmb, ...
                       mcData.inputs.At, mcData.inputs.Ae, mcData.inputs.dExt, ...
                       mcData.inputs.dInt, mcData.inputs.L, mcData.inputs.alpha];
        inputNames = {'Burn Rate Coeff (a)', 'Burn Rate Exp (n)', 'O/F Ratio', 'Ambient Pressure', ...
                    'Throat Area (At)', 'Exit Area (Ae)', 'Grain OD (dExt)', ...
                    'Port Dia (dInt)', 'Grain Length (L)', 'Div Angle (\alpha)'};

        if config.analysisMode == 0 || config.analysisMode == 2
            % MEOP Tornado Plot
            corrVals = zeros(length(inputNames), 1);
            for i = 1:length(inputNames)
                if std(inputMatrix(:,i)) == 0 % Properly catches frozen variables
                    corrVals(i) = 0;
                else
                    corrVals(i) = corr(inputMatrix(:,i), mcData.MEOP, 'Type', 'Spearman');
                end
            end
    
            [~, sortIdx] = sort(abs(corrVals), 'ascend');
            sortedNames = inputNames(sortIdx); actualVals = corrVals(sortIdx);
    
            figure('Color', 'w', 'Name', 'Tornado Plot - MEOP Sensitivity Drivers');
            b = barh(actualVals, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.7);
            for k = 1:length(actualVals)
                if actualVals(k) > 0; b.CData(k,:) = techOrange; else; b.CData(k,:) = techPurple; end
            end
            yticks(1:length(sortedNames)); yticklabels(sortedNames);
            xlabel('Spearman Correlation Coefficient [-]', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
            title(sprintf('MEOP Drivers (Mode %d)', analysisMode), 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
            grid on; set(gca, 'FontName', fontName, 'FontSize', 14, 'GridLineStyle', '--', 'GridAlpha', 0.3);
            
            
    
            % Burn Time Tornado Plot
            corrValsBT = zeros(length(inputNames), 1);
            for i = 1:length(inputNames)
                if std(inputMatrix(:,i)) == 0 % Catch frozen variables
                    corrValsBT(i) = 0;
                else
                    corrValsBT(i) = corr(inputMatrix(:,i), mcData.burnTime, 'Type', 'Spearman');
                end
            end
    
            [~, sortIdxBT] = sort(abs(corrValsBT), 'ascend');
            sortedNamesBT = inputNames(sortIdxBT); 
            actualValsBT = corrValsBT(sortIdxBT);
    
            figure('Color', 'w', 'Name', 'Tornado Plot - Burn Time Sensitivity Drivers');
            bBT = barh(actualValsBT, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.7);
            for k = 1:length(actualValsBT)
                if actualValsBT(k) > 0
                    bBT.CData(k,:) = techOrange; 
                else
                    bBT.CData(k,:) = techPurple;  
                end
            end
            yticks(1:length(sortedNamesBT)); yticklabels(sortedNamesBT);
            xlabel('Spearman Correlation Coefficient [-]', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
            title(sprintf('What Drives Burn Time? (Mode %d)', analysisMode), 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
            grid on; set(gca, 'FontName', fontName, 'FontSize', 14, 'GridLineStyle', '--', 'GridAlpha', 0.3);
        end

        % Trade-off Matrix
        Y = [mcData.MEOP(:), mcData.TMax(:), mcData.IspAvg(:), mcData.burnTime(:)];
        Rmat = corrcoef(Y, 'Rows', 'complete');
        figure('Color', 'w', 'Name','Monte Carlo - Output Correlation');
        imagesc(Rmat); colorbar; axis equal tight;
        set(gca, 'XTick', 1:4, 'XTickLabel', {'MEOP','Tmax','Isp_{avg}','Burn Time'}, 'FontName', fontName)
        set(gca, 'YTick', 1:4, 'YTickLabel', {'MEOP','Tmax','Isp_{avg}','Burn Time'}, 'FontName', fontName)
        title('Trade-Off Matrix (Output Coupling)', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold')

        % Joint Envelope Trade-Off
        figure('Color', 'w', 'Name', 'Joint Envelope: MEOP vs Max Thrust');
        hold on; box on; grid on;
        scatter(mcData.MEOP * 1e-5, mcData.TMax * 1e-3, 25, techBlue, 'filled', 'MarkerFaceAlpha', 0.5);
        xline(mcData.MEOPMean * 1e-5, '--', 'Mean MEOP', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
        yline(mean(mcData.TMax) * 1e-3, '--', 'Mean Thrust', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
        xlabel('MEOP [bar]', 'FontWeight', 'bold'); ylabel('Maximum Thrust [kN]', 'FontWeight', 'bold');
        title('System Design Envelope (MEOP vs Thrust)', 'FontWeight', 'bold');
        set(gca, 'FontName', fontName, 'FontSize', fontSize, 'GridLineStyle', '--', 'GridAlpha', 0.3);
    end

    %% Sensitivities (OAT & Mode 0/1 Bar Charts)
    if analysisMode == 0 || analysisMode == 1
        
        % Deterministic Sensitivity Magnitude (OAT Deltas)
        % Calculate the absolute peak-to-peak spread caused by the +/- 3 sigma sweep
        deltaMeopA = (max(mcData.oat.MEOPVaryA) - min(mcData.oat.MEOPVaryA)) * 1e-5; % in bar
        deltaMeopN = (max(mcData.oat.MEOPVaryN) - min(mcData.oat.MEOPVaryN)) * 1e-5;
        
        deltaBtA = max(mcData.oat.BTVaryA) - min(mcData.oat.BTVaryA); % in seconds
        deltaBtN = max(mcData.oat.BTVaryN) - min(mcData.oat.BTVaryN);

        figure('Color', 'w', 'Name', 'Simplified Sensitivity: a & n Drivers');
        
        % Because MEOP (bar) and Burn Time (s) are completely different units, 
        % is used a tiled layout or subplots to show their pure physical deltas.
        subplot(1, 2, 1); hold on; grid on; box on;
        b1 = bar([deltaMeopA, deltaMeopN], 'FaceColor', 'flat', 'EdgeColor', 'none');
        b1.CData(1,:) = techBlue; b1.CData(2,:) = techOrange;
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Coefficient (a)', 'Exponent (n)'}, 'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('\Delta MEOP Spread [bar]', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
        title('Max Pressure Impact', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
        set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);

        subplot(1, 2, 2); hold on; grid on; box on;
        b2 = bar([deltaBtA, deltaBtN], 'FaceColor', 'flat', 'EdgeColor', 'none');
        b2.CData(1,:) = techBlue; b2.CData(2,:) = techOrange;
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Coefficient (a)', 'Exponent (n)'}, 'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('\Delta Burn Time Spread [s]', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
        title('Max Burn Time Impact', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
        set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
        
        sgtitle('Absolute Deterministic Impact (\pm 3\sigma Sweep)', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');

        % Marginal Sensitivity Scatter Plot
        figure('Color', 'w', 'Name', 'OAT Marginal Sensitivity: a vs n');
        
        % MEOP OAT
        subplot(1, 2, 1); hold on; grid on; box on;
        plot(mcData.oat.pctDevA, mcData.oat.MEOPVaryA * 1e-5, '-', 'Color', techBlue, 'LineWidth', 2.5, 'DisplayName', 'Varying only ''a''');
        plot(mcData.oat.pctDevN, mcData.oat.MEOPVaryN * 1e-5, '-', 'Color', techOrange, 'LineWidth', 2.5, 'DisplayName', 'Varying only ''n''');
        plot(0, mcData.nominals.simMEOP * 1e-5, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Simulated Nominal');
        xlabel('Parameter Deviation [%]', 'FontWeight', 'bold'); ylabel('MEOP [bar]', 'FontWeight', 'bold');
        title('Isolated Impact on MEOP', 'FontWeight', 'bold');
        legend('Location', 'northwest'); 
        set(gca, 'FontName', fontName, 'FontSize', 14, 'GridAlpha', 0.3, 'GridLineStyle', '--');
    
        % Burn Time OAT
        subplot(1, 2, 2); hold on; grid on; box on;
        plot(mcData.oat.pctDevA, mcData.oat.BTVaryA, '-', 'Color', techBlue, 'LineWidth', 2.5, 'DisplayName', 'Varying only ''a''');
        plot(mcData.oat.pctDevN, mcData.oat.BTVaryN, '-', 'Color', techOrange, 'LineWidth', 2.5, 'DisplayName', 'Varying only ''n''');
        plot(0, mcData.nominals.simBurnTime, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Simulated Nominal');
        xlabel('Parameter Deviation [%]', 'FontWeight', 'bold'); ylabel('Burn Time [s]', 'FontWeight', 'bold');
        title('Isolated Impact on Burn Time', 'FontWeight', 'bold');
        legend('Location', 'northeast'); 
        set(gca, 'FontName', fontName, 'FontSize', 14, 'GridAlpha', 0.3, 'GridLineStyle', '--');
    end

    %% Cooling Jacket Risk
    if config.runCoolingJacket && analysisMode == 1
        % Extract Data
        TOutlet = mcData.thermal.TOutlet;
        maxQ = mcData.thermal.maxQ;
        boilingLimit = mcData.thermal.boilingLimit;
        
        TBCHot = mcData.thermal.TBCHot;
        TBCMaxLimit = mcData.thermal.TBCMaxLimit;
        
        INCONELHot = mcData.thermal.INCONELHot;
        INCONELMaxLimit = mcData.thermal.INCONELMaxLimit;
        
        % Bundle plot config for the helper function
        pConfig.fontName = fontName;
        pConfig.fontSize = fontSize;
        pConfig.techOrange = techOrange;
        
        figure('Color', 'w', 'Name', 'Cooling Jacket Risk Analysis', 'WindowState', 'maximized');
        sgtitle(sprintf('Cooling Jacket Thermal Analysis for water mass flow rate %.2f kg/s', mcData.thermal.mDotWater), ...
            'FontName', fontName, 'FontSize', fontSize + 2, 'FontWeight', 'bold');

        % Subplot 1: TOutlet vs Boiling Limit
        subplot(2, 2, 1);
        plotThermalRisk(TOutlet, boilingLimit, 'Boiling Limit', 'Max Coolant Wall Temp [K]', 'Coolant Thermal Risk', techBlue, pConfig);

        % Subplot 2: maxQ Distribution (No explicit limit line)
        subplot(2, 2, 2);
        hold on; grid on; box on;
        histogram(maxQ, 'Normalization', 'pdf', 'FaceColor', techOrange, 'EdgeColor', 'w');
        plotNormalFit(maxQ);
        dataSpreadQ = max(maxQ) - min(maxQ);
        if dataSpreadQ > 0
            xlim([min(maxQ) - (dataSpreadQ*0.5), max(maxQ) + (dataSpreadQ*0.5)]); 
        end
        xline(mean(maxQ), '--', 'Mean Heat Flux', 'Color', [0.3 0.3 0.3], 'LineWidth', 2, 'LabelVerticalAlignment', 'top', 'FontName', fontName, 'FontSize', 14);
        xlabel('Maximum Heat Flux (maxQ) [W/m^2]', 'FontWeight', 'bold');
        ylabel('Probability Density [-]', 'FontWeight', 'bold');
        title('Cooling Jacket Heat Flux Distribution', 'FontWeight', 'bold');
        set(gca, 'FontName', fontName, 'FontSize', fontSize, 'GridLineStyle', '--');
        
        % Subplot 3: TBC Temp vs Limit
        subplot(2, 2, 3);
        plotThermalRisk(TBCHot, TBCMaxLimit, 'TBC Limit', 'TBC Maximum Temperature [K]', 'TBC Thermal Risk', techPurple, pConfig);

        % Subplot 4: INCONEL Temp vs Limit
        subplot(2, 2, 4);
        plotThermalRisk(INCONELHot, INCONELMaxLimit, 'INCONEL Limit', 'INCONEL Maximum Temperature [K]', 'INCONEL Thermal Risk', techGreen, pConfig);
    end
    % pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom pom 
    %% Local Helper Functions 
    
    function plotRobustHistogram(data, faceColor, titleText, fontName)
        % plotRobustHistogram Handles plotting gracefully when variance is near zero.
        mu = mean(data, 'omitnan');
        sig = std(data, 'omitnan');
        
        % Check for zero or near-zero variance (relative and absolute checks)
        if sig < 1e-6 * abs(mu) || sig < 1e-9
            % Plot as a probability block to avoid PDF scaling singularities
            histogram(data, 1, 'Normalization', 'probability', 'FaceColor', faceColor, 'EdgeColor', 'k');
            
            % Force breathing room on the X-axis so it isn't completely flush
            if mu == 0
                xlim([-1, 1]);
            else
                xlim([mu * 0.9, mu * 1.1]);
            end
            title(titleText, 'FontName', fontName);
        else
            % Standard PDF plot with Gaussian fit overlay
            histogram(data, 30, 'Normalization', 'pdf', 'FaceColor', faceColor);
            plotNormalFit(data);
            title(titleText, 'FontName', fontName);
        end
    end

    function plotNormalFit(data)
        mu = mean(data, 'omitnan');
        sig = std(data, 'omitnan');
        xRange = linspace(min(data), max(data), 200);
        yFit = (1 / (sig * sqrt(2*pi))) * exp(-0.5 * ((xRange - mu) / sig).^2);
        hold on;
        plot(xRange, yFit, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Normal Fit');
    end
    
    function plotThermalRisk(data, limitVal, limitName, xLabelText, titleText, faceColor, pConfig)
        % plotThermalRisk Handles the repetitive histogram + dynamic limit line logic
        hold on; grid on; box on;
        histogram(data, 'Normalization', 'pdf', 'FaceColor', faceColor, 'EdgeColor', 'w');
        plotNormalFit(data);
        
        % Bounds and formatting
        dataSpread = max(data) - min(data);
        if dataSpread > 0
            xlim([min(data) - (dataSpread*0.5), max(data) + (dataSpread*0.5)]); 
        end
        xline(mean(data), '--', 'Mean', 'Color', [0.3 0.3 0.3], 'LineWidth', 2, 'LabelVerticalAlignment', 'top', 'FontName', pConfig.fontName, 'FontSize', 14);
        
        % Dynamic Limit Line logic
        xl = xlim;
        yl = ylim;
        if limitVal >= xl(1) && limitVal <= xl(2)
            xline(limitVal, '-.', limitName, 'Color', pConfig.techOrange, 'LineWidth', 2.5, ...
                'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left', ...
                'FontName', pConfig.fontName, 'FontSize', 14, 'FontWeight', 'bold');
        elseif limitVal > xl(2)
            text(xl(2), yl(2)*0.85, sprintf('%s (%.1f K) \\rightarrow  ', limitName, limitVal), ...
                'Color', pConfig.techOrange, 'FontWeight', 'bold', 'FontSize', 14, 'FontName', pConfig.fontName, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        else
            text(xl(1), yl(2)*0.85, sprintf('  \\leftarrow %s (%.1f K)', limitName, limitVal), ...
                'Color', pConfig.techOrange, 'FontWeight', 'bold', 'FontSize', 14, 'FontName', pConfig.fontName, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end
        
        % Labels and Probability calculation
        xlabel(xLabelText, 'FontWeight', 'bold');
        ylabel('Probability Density [-]', 'FontWeight', 'bold');
        title(titleText, 'FontWeight', 'bold');
        set(gca, 'FontName', pConfig.fontName, 'FontSize', pConfig.fontSize, 'GridLineStyle', '--');
        
        % Localized Probability Box
        PoF = (sum(data >= limitVal) / length(data)) * 100;
        text(xl(1) + (xl(2)-xl(1))*0.03, yl(2)*0.95, sprintf('P(Exceed): %.2f%%', PoF), ...
            'BackgroundColor', 'w', 'EdgeColor', pConfig.techOrange, 'FontName', pConfig.fontName, ...
            'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    end
end