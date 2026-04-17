% [name] - this script provides nominal performances computed through a
% Monte Carlo simulation to consider uncertainties on the Vieille's law
% parameters' uncertainty

%   Available parameters:
%      - a: Vieille's law coefficient
%      - n: Vieille's law exponent

%% Setup

clear; close all; clc

addpath(genpath("./src"))
plotFlag = 1;

if ~exist(fullfile('.', 'data'), 'dir') || ~exist(fullfile('.', 'data', 'propellant.mat'), 'file')
    data = savePropellantData();
else
    data = load(fullfile('.', 'data', 'propellant.mat'));
end

% calling all the variables/structures coming from the main
main;



iter = 1000;
uncFuncChoice = 2;

aNom = a;
nNom = n;

%%
switch uncFuncChoice
    case 1
        % --- CASE 1: INDEPENDENT SAMPLING WITH MESHGRID ---
        [~, aSigmaMC, ~, nSigma, ~] = uncertaintyVieille(data.ccPressure * 1e-6, data.burnRate * 1e3);
        aSigmaMC = aSigmaMC * 1e-3/(10^(6*nNom));
        
        aMCunshuffled = aNom * ones(iter, 1) + randn(iter, 1) * aSigmaMC;
        nMCunshuffled = nNom * ones(iter, 1) + randn(iter, 1) * nSigma;
        
        % initializing the shuffle => incorporated with meshgrid
        [aMC_matrix, nMC_matrix] = meshgrid(aMCunshuffled, nMCunshuffled);
        
        % flatten the matrices into 1D vectors for the parfor loop
        aMC = aMC_matrix(:);
        nMC = nMC_matrix(:);
        
    case 2
        % --- CASE 2: CORRELATED MULTIVARIATE SAMPLING ---
        [~, ~, ~, mu_qm, Sigma_qm] = uncertaintyVieilleDos(data.ccPressure * 1e-6, data.burnRate * 1e3);
        
        % To match the size of the meshgrid output (iter * iter), we directly 
        % draw that many samples from the correlated multivariate distribution.
        total_samples = iter^2; 
        samples_qm = mvnrnd(mu_qm, Sigma_qm, total_samples);
        
        qMC = samples_qm(:, 1);
        nMC = samples_qm(:, 2);
        
        % Transform q back to a and convert to SI units
        aMC = exp(qMC);
        aMC = aMC .* 1e-3 ./ (10.^(6 .* nMC)); 
        
    otherwise
        error("flag not defined.")
end


%%


rhoP = data.cea.rhoP;
cStar = performanceNom.cstar;
Athroat = nozzle.At;

% Use numel to get the total number of combinations (100 * 100 = 10000)
iterMC = numel(aMC);

burntime = zeros(iterMC, 1);
MEOP = zeros(iterMC, 1);

%%

parfor index = 1:iterMC
    a = aMC(index);
    n = nMC(index);
    [t, p, rb] = computeBurn(a, n, rhoP, cStar, grain, Athroat);
    burntime(index) = max(t);
    MEOP(index) = max(p);
end

%%

meanBT = mean(burntime);
stdBT = std(burntime);
meanMEOP = mean(MEOP);
stdMEOP = std(MEOP);


% Cumulative tracking arrays
progBTMean = cumsum(burntime) ./ (1:iterMC)';
progBTDevStd = sqrt((cumsum(burntime.^2) ./ (1:iterMC)') - progBTMean.^2);

progMEOPMean = cumsum(MEOP) ./ (1:iterMC)';
progMEOPDevStd = sqrt((cumsum(MEOP.^2) ./ (1:iterMC)') - progMEOPMean.^2);

% Relative Error Calculations (%) compared to the final aggregated value
errorBTMean = abs(progBTMean - meanBT) / meanBT * 100;
errorBTDevStd = abs(progBTDevStd - stdBT) / stdBT * 100;

errorMEOPMean = abs(progMEOPMean - meanMEOP) / meanMEOP * 100;
errorMEOPDevStd = abs(progMEOPDevStd - stdMEOP) / stdMEOP * 100;

%% Plotting Relative Convergence


figure('Name', 'Relative Error Convergence Monitoring');

% Burn Time Mean Error
subplot(2, 2, 1);
semilogy(1:iterMC-1, errorBTMean(1:end-1), 'b-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Relative Error (%)');
title('Burn Time: Mean Convergence');
grid on;

% Burn Time Standard Deviation Error
subplot(2, 2, 2);
semilogy(1:iterMC-1, errorBTDevStd(1:end-1), 'r-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Relative Error (%)');
title('Burn Time: Std Dev Convergence');
grid on;

% MEOP Mean Error
subplot(2, 2, 3);
semilogy(1:iterMC-1, errorMEOPMean(1:end-1), 'b-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Relative Error (%)');
title('MEOP: Mean Convergence');
grid on;

% MEOP Standard Deviation Error
subplot(2, 2, 4);
semilogy(1:iterMC-1, errorMEOPDevStd(1:end-1), 'r-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Relative Error (%)');
title('MEOP: Std Dev Convergence');
grid on;


%% Plotting

figure('Name', 'Convergence Monitoring');

% Burn Time Mean
subplot(2, 2, 1);
semilogy(1:iterMC, progBTMean, 'b-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Burn Time (s)');
title('Mean Burn Time');
grid on;

% Burn Time Standard Deviation
subplot(2, 2, 2);
plot(1:iterMC, progBTDevStd, 'r-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('\sigma Burn Time (s)');
title('Standard Deviation of Burn Time');
grid on;

% MEOP Mean
subplot(2, 2, 3);
plot(1:iterMC, progMEOPMean, 'b-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('MEOP (Pa)');
title('Mean MEOP');
grid on;

% MEOP Standard Deviation
subplot(2, 2, 4);
plot(1:iterMC, progMEOPDevStd, 'r-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('\sigma MEOP (Pa)');
title('Standard Deviation of MEOP');
grid on;