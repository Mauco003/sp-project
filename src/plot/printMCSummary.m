function printMCSummary(mcData, config)
% printMCSummary -  Generates a clean, professional command window printout
% summarizing the key statistical and deterministic results of the simulation.
%
% SYNTAX:
%  printMCSummary(mcData, config)
%
% INPUT:
%  mcData - A structure containing all relevant Monte Carlo simulation results and statistics.
%  config - A structure containing configuration settings for the analysis

    % Helper for percentage calculations (Coefficient of Variation)
    covMEOP = (mcData.MEOPStd / mcData.MEOPMean) * 100;
    covBT   = (mcData.burnTimeStd / mcData.burnTimeMean) * 100;
    
    fprintf('\n');
    fprintf('===========================================================================\n');
    fprintf('                 MONTE CARLO SIMULATION SUMMARY\n');
    fprintf('===========================================================================\n');
    fprintf(' Total Runs       : %d (Mode %d)\n', mcData.totalRunsMC, config.analysisMode);
    
    if config.analysisMode == 0 || config.analysisMode == 1
        fprintf(' Sampling Type    : %s\n', upper(config.correlationType));
    end
    fprintf('---------------------------------------------------------------------------\n');
    
    %% Primary Output Statistics
    fprintf(' [PRIMARY BALLISTICS]\n');
    fprintf('                    Sim. Nominal |      Mean      |     3σ Range    \n');
    
    % MEOP Row (Converted to bar)
    nomMeop = mcData.nominals.simMEOP * 1e-5;
    muMeop  = mcData.MEOPMean * 1e-5;
    sigMeop = mcData.MEOPStd * 1e-5;
    fprintf(' MEOP (bar)       : %12.2f | %6.2f (±%4.1f%%) | [%.1f, %.1f]\n', ...
        nomMeop, muMeop, covMEOP, muMeop - 3*sigMeop, muMeop + 3*sigMeop);
    
    % Burn Time Row
    nomBt = mcData.nominals.simBurnTime;
    muBt  = mcData.burnTimeMean;
    sigBt = mcData.burnTimeStd;
    fprintf(' Burn Time (s)    : %12.2f | %6.2f (±%4.1f%%) | [%.1f, %.1f]\n', ...
        nomBt, muBt, covBT, muBt - 3*sigBt, muBt + 3*sigBt);
    fprintf('---------------------------------------------------------------------------\n');
    
    %% Comprehensive Performance
    if config.useComprehensive
        fprintf(' [COMPREHENSIVE PERFORMANCE]\n');
        
        % Convert Thrust to kN
        muTmax = mean(mcData.TMax) * 1e-3;
        peakTmax = max(mcData.TMax) * 1e-3;
        muTavg = mean(mcData.TAvg) * 1e-3;
        
        fprintf(' Max Thrust (kN)  : Mean = %7.2f | Absolute Peak = %7.2f\n', muTmax, peakTmax);
        fprintf(' Avg Thrust (kN)  : Mean = %7.2f \n', muTavg);
        
        % Note: Isp depends on your units (usually m/s or s). Printed raw.
        fprintf(' Isp Avg          : Mean = %7.2f \n', mean(mcData.IspAvg));
        
        % Reliability
        avgRel = mean(mcData.thrustReliability, 'omitnan');
        fprintf(' Reliability      : On average, thrust is >90%% nominal for %.1f%% of burn.\n', avgRel);
        
        fprintf('---------------------------------------------------------------------------\n');
    end

    %% Deterministic OAT Sensitivity
    if config.analysisMode == 0 || config.analysisMode == 1
        % Calculate absolute Peak-to-Peak Delta from the deterministic sweep
        deltaMeopA = (max(mcData.oat.MEOPVaryA) - min(mcData.oat.MEOPVaryA)) * 1e-5;
        deltaMeopN = (max(mcData.oat.MEOPVaryN) - min(mcData.oat.MEOPVaryN)) * 1e-5;
        
        deltaBtA = max(mcData.oat.BTVaryA) - min(mcData.oat.BTVaryA);
        deltaBtN = max(mcData.oat.BTVaryN) - min(mcData.oat.BTVaryN);

        fprintf(' [DETERMINISTIC SENSITIVITY (± 3σ Sweep Peak-to-Peak Deltas)]\n');
        fprintf('                    Varying ''a''   |   Varying ''n'' \n');
        fprintf(' Δ MEOP (bar)     : %12.2f | %13.2f \n', deltaMeopA, deltaMeopN);
        fprintf(' Δ Burn Time (s)  : %12.2f | %13.2f \n', deltaBtA, deltaBtN);
        fprintf('---------------------------------------------------------------------------\n');
    end

    %% Thermal Risk
    if config.runCoolingJacket && config.analysisMode == 1
        fprintf(' [COOLING JACKET RISK ANALYSIS (water mass flow rate %f kg/s)]\n', mcData.thermal.mDotWater);
        muTemp = mean(mcData.thermal.TOutlet, 'omitnan');
        maxTemp = max(mcData.thermal.TOutlet, [], 'omitnan');
        
        % Calculate maxQ stats
        muMaxQ = mean(mcData.thermal.maxQ, 'omitnan');
        peakMaxQ = max(mcData.thermal.maxQ, [], 'omitnan');
        
        fprintf(' Maximum Coolant Side Wall temperature (K)  : Mean = %6.2f | Max Observed  = %6.2f\n', muTemp, maxTemp);
        fprintf(' Max Heat Flux    : Mean = %10.2e | Max Observed  = %10.2e\n', muMaxQ, peakMaxQ);
        fprintf(' Boiling Limit    : %6.2f K\n', mcData.thermal.boilingLimit);
        
        if mcData.thermal.PoF > 0
            fprintf(' STATUS           : WARNING - Probability of Failure = %.2f%%\n', mcData.thermal.PoF);
        else
            fprintf(' STATUS           : SAFE - 0%% probability of localized boiling.\n');
        end
        fprintf('===========================================================================\n\n');
    else
        fprintf('===========================================================================\n\n');
    end
end

