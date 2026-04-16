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


[~, aSigma, ~, nSigma, ~] = uncertaintyVieille(data.ccPressure * 1e-5, data.burnRate * 1e3);

aNom = a;
nNom = n;

aSigma = aSigma * 1e-3/(10^(5*n));

rhoP = data.cea.rhoP;
cStar = performanceNom.cstar;
Athroat = nozzle.At;

iter = 100;
aMCunshuffled = aNom * ones(iter, 1) + randn(iter, 1) * aSigma;
nMCunshuffled = nNom * ones(iter, 1) + randn(iter, 1) * nSigma;

% initializing the shuffle => incorporated with meshgrid

[aMC, nMC] = meshgrid(aMCunshuffled, nMCunshuffled);

matrixMC = [aMC(:), nMC(:)];
iterMC = size(matrixMC, 1);

burntime = zeros(iterMC, 1);
MEOP = zeros(iterMC, 1);

for index = 1:iterMC

    a = matrixMC(index, 1);
    n = matrixMC(index, 2);
    [t, p, rb] = computeBurn(a, n, rhoP, cStar, grain, Athroat);
    burntime(index) = max(t);
    MEOP(index) = max(p);

end

meanBT = mean(burntime);
stdBT = std(burntime);
meanMEOP = mean(MEOP);
stdMEOP = std(MEOP);

progBTMean = cumsum(burntime) ./ (1:iterMC)';
progBTDevStd = sqrt((cumsum(burntime.^2) ./ (1:iterMC)') - progBTMean.^2);

errorBTMean = abs(meanBT * ones(iterMC, 1) - progBTMean)/meanBT * 100;

figure()
semilogy(1:iterMC-1, errorBTMean(1:end-1))

errorBTDevStd = abs(stdBT * ones(iterMC, 1) - progBTDevStd)/stdBT * 100;

figure()
semilogy(1:iterMC-1, errorBTDevStd(1:end-1))


progMEOPMean = cumsum(MEOP) ./ (1:iterMC)';
progMEOPDevStd = sqrt((cumsum(MEOP.^2) ./ (1:iterMC)') - progMEOPMean.^2);


% Plotting the results
figure;
subplot(2, 1, 1);
plot(1:iterMC, progBTMean, 'b-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Burn Time (s)');
title('Mean Burn Time with Standard Deviation');
grid on;