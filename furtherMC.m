% furtherMC - This script runs a Monte Carlo simulation to 
%   analyze the performance of a rocket engine design under 
%   varying conditions. 
%
%   Available parameters:
%      - a: Vieille's law coefficient
%      - n: Vieille's law exponent
%      - O/F: Oxidizer to fuel mass ratio
%      - pa: Ambient pressure
%      - At: Nozzle throat area
%      - Ae: Nozzle exit area
%      - dInt0: Initial grain internal diameter
%      - dExt0: Initial grain external diameter
%      - L0: Initial grain length
%      - mPTot: Total propellant mass
%      - alpha: Nozzle divergent half-angle
%      - beta: Nozzle convergent half-angle

%%

main;

%% standard deviations
% a & n already known

% OF mass ratio
OFNominal = 80/20;
% excursion between +-2% 
OFSigma = (1.02*OFNominal - OFNominal)/3; % standard deviation for O/F ratio

% ambient pressure
% excursions between 97000 Pa e 105000 Pa (to verify) => +-4% variation
pambSigma = (1.04*constants.pAmb - constants.pAmb)/3;

% throat diameter => uncertainty on the diameter
% excursion between +-1%
dtNominal = sqrt(4/pi * nozzle.At);
dtSigma = (1.01 * dtNominal - dtNominal)/3;

% exit area => uncertainty on the diameter
% excursion between +-1%
deNominal = sqrt(4/pi * nozzle.Ae);
deSigma = (1.01 * deNominal - deNominal)/3;

% external and internal diameter of the grain
% excursion between +-0.5%
dextSigma = (1.005 * grain.dExt0  - grain.dExt0)/3;

dintSigma = (1.005 * grain.dInt0  - grain.dInt0)/3;

% initial length of the grain
% same excursion as the diameters I guess
LSigma = (1.005 * grain.L0 - grain.L0)/3;

% propellant mass
% we can evaluate it as a combination of an uncertainty on the density of
% the propellant and an uncertainty on the volume
% for now I'll set +-2% 
mPSigma = (1.02 * mPTot - mPTot)/3;

% alpha & beta
L_conv = nozzle.x(2);
beta = 30; % deg
L_parete_conv = L_conv/cosd(beta);
delta_parete = 0.00001;
delta_beta = atand(delta_parete/L_parete_conv);
betaSigma  = delta_beta/3;

L_div = nozzle.x(3);
alpha = 15;
L_parete_div = L_div/cosd(alpha);
delta_alpha = atand(delta_parete/L_parete_div);
alphaSigma = delta_alpha/3;

%% Initializing the Monte Carlo simulation

iter = 100;
uncFuncChoice = 2;


aNominal = a;
nNominal = n;

switch uncFuncChoice
    case 1
        % Independent sampling with meshgrid
        [~, aSigmaMC, ~, nSigma, ~] = uncertaintyVieille(propellant.ccPressure * 1e-6, propellant.burnRate * 1e3);
        aSigmaMC = aSigmaMC * 1e-3/(10^(6*nNominal));
        
        aMCunshuffled = aNominal * ones(iter, 1) + randn(iter, 1) * aSigmaMC;
        nMCunshuffled = nNominal * ones(iter, 1) + randn(iter, 1) * nSigma;
        
        % Initializing the shuffle => incorporated with meshgrid
        [aMC_matrix, nMC_matrix] = meshgrid(aMCunshuffled, nMCunshuffled);
        
        % Flatten the matrices into 1D vectors
        aMC_array = aMC_matrix(:);
        nMC_array = nMC_matrix(:);
        
    case 2
        % Correlated multivariate sampling
        % The physics behind the problem teaches us that is nearly impossible having 
        % propellant with combination of a and n s.t. the performace have maximum values 
        % much much higher. The regression has a way of expressing the relation between 
        % the two with the covariance. So this should be exploited ass well
        [~, ~, ~, mu_qm, Sigma_qm] = uncertaintyVieilleDos(propellant.ccPressure * 1e-6, propellant.burnRate * 1e3);
        
        % To match the size of the meshgrid output (iter * iter)
        total_samples = iter^2; 
        samples_qm = mvnrnd(mu_qm, Sigma_qm, total_samples);
        
        qMC = samples_qm(:, 1);
        nMC_array = samples_qm(:, 2);
        
        % Transform q back to a and convert to SI units
        aMC_array = exp(qMC);
        aMC_array = aMC_array .* 1e-3 ./ (10.^(6 .* nMC_array)); 
        
    otherwise
        error("flag not defined.")
end

mcmax = numel(aMC_array); % Total number of combinations

propellantMC = propellant;
constantsMC = constants;
nozzleMC = nozzle;
grainMC = grain;


% Pre-generate all random inputs outside the loop (vectorized)
aVec    = reshape(aMC_array, 1, mcmax);
nVec    = reshape(nMC_array, 1, mcmax);
OFVec   = OFNominal + randn(1, mcmax) * OFSigma;
pambVec = constants.pAmb + randn(1, mcmax) * pambSigma;
AtVec   = (dtNominal + randn(1, mcmax) * dtSigma).^2 * pi/4;
AeVec   = (deNominal + randn(1, mcmax) * deSigma).^2 * pi/4;
betaVec = beta + randn(1, mcmax) * betaSigma;
alphaVec= alpha + randn(1, mcmax) * alphaSigma;
dextVec = grain.dExt0 + randn(1, mcmax) * dextSigma;
dintVec = grain.dInt0 + randn(1, mcmax) * dintSigma;
LVec    = grainMC.L0 + randn(1, mcmax) * LSigma;
MpVec   = mPTot + randn(1, mcmax) * mPSigma;


% Pre-allocate output arrays
T_max = zeros(1,mcmax);
T_avg = zeros(1,mcmax);
Isp_max = zeros(1,mcmax);
Isp_avg = zeros(1,mcmax);
time_pc = zeros(1,mcmax);

rB = zeros(1,mcmax);
Isp_tot_max = zeros(1,mcmax);
deltatGuaranteedThrust = zeros(1,mcmax);

%% Physics Execution Loop

for index = 1:mcmax
    % Extract the pre-calculated random values for this iteration
    aMc    = aVec(index);
    nMc    = nVec(index);
    OFMc   = OFVec(index);
    pambMc = pambVec(index);
    AtMc   = AtVec(index);
    AeMc   = AeVec(index);
    betaMc = betaVec(index); 
    alphaMc= alphaVec(index);
    dextMc = dextVec(index);
    dintMc = dintVec(index);
    LMc    = LVec(index);
    MpMc   = MpVec(index);

    % Setup the structures for this run
    rB(index) = aMc * input.pcNominal^nMc;

    MOx = OFMc/(1 + OFMc);
    MF  = 1/(1 + OFMc);
    propellantMC.cea.rhoP = 1/(MOx/propellant.cea.rhoAP + MF/propellant.cea.rhoHTPB);
    
    [ChamberData, ~] = getThermoProfileCEA_froz(MOx*1e2, MF*1e2, input.pcNominal, 2, 2);
    propellantMC.cea.gamma = ChamberData.gamma;

    constantsMC.pAmb = pambMc;
    nozzleMC.At = AtMc;
    nozzleMC.Ae = AeMc;

    grainMC.dExt0 = dextMc;
    grainMC.dInt0 = dintMc;
    grainMC.L0    = LMc;

    % Run computations
    [tMC_output, performanceMC_output, grainMC_output] = ...
        computePerformance(aMc, nMc, propellantMC, performanceNom, nozzleMC, grainMC, constantsMC);

    I_tot = trapz(tMC_output, performanceMC_output.thrust);

    % Store outputs
    T_max(index) = max(performanceMC_output.thrust);
    T_avg(index) = mean(performanceMC_output.thrust);
    Isp_avg(index) = mean(performanceMC_output.Isp);
    Isp_max(index) = max(performanceMC_output.Isp);
    Isp_tot_max(index) = I_tot;

    threshold = 0.9 * max(performanceMC_output.thrust);
    idx = performanceMC_output.thrust >= threshold;

    if any(idx)
        timevector = tMC_output(idx);
        deltatGuaranteedThrust(index) = timevector(end) - timevector(1);
    end
end

%% Monte Carlo post-processing plots

% Summary statistics
fprintf('\nMonte Carlo summary (%d runs)\n', mcmax);
fprintf('T_max   : mean = %.3f, std = %.3f\n', mean(T_max), std(T_max));
fprintf('T_avg   : mean = %.3f, std = %.3f\n', mean(T_avg), std(T_avg));
fprintf('Isp_max : mean = %.3f, std = %.3f\n', mean(Isp_max), std(Isp_max));
fprintf('Isp_avg : mean = %.3f, std = %.3f\n', mean(Isp_avg), std(Isp_avg));
fprintf('Impulse : mean = %.3f, std = %.3f\n', mean(Isp_tot_max), std(Isp_tot_max));
fprintf('t_90%%   : mean = %.3f, std = %.3f\n', mean(deltatGuaranteedThrust), std(deltatGuaranteedThrust));

% Histograms of main outputs
figure('Name','Monte Carlo - Main Output Distributions');

subplot(2,3,1)
histogram(T_max, 30)
xlabel('Peak thrust')
ylabel('Count')
title('T_{max}')

subplot(2,3,2)
histogram(T_avg, 30)
xlabel('Average thrust')
ylabel('Count')
title('T_{avg}')

subplot(2,3,3)
histogram(Isp_max, 30)
xlabel('Peak Isp')
ylabel('Count')
title('Isp_{max}')

subplot(2,3,4)
histogram(Isp_avg, 30)
xlabel('Average Isp')
ylabel('Count')
title('Isp_{avg}')

subplot(2,3,5)
histogram(Isp_tot_max, 30)
xlabel('Total impulse')
ylabel('Count')
title('Total impulse')

subplot(2,3,6)
histogram(deltatGuaranteedThrust, 30)
xlabel('Time above 90% max thrust [s]')
ylabel('Count')
title('\Delta t guaranteed thrust')

sgtitle('Monte Carlo Output Distributions')

% Burn-rate distribution
figure('Name','Monte Carlo - Burn Rate');
histogram(rB, 30)
xlabel('Burn rate')
ylabel('Count')
title('Distribution of burn rate r_B')

% Boxplots for quick comparison
figure('Name','Monte Carlo - Boxplots');
boxplot([T_max(:), T_avg(:), Isp_max(:), Isp_avg(:), Isp_tot_max(:), deltatGuaranteedThrust(:)], ...
    'Labels', {'Tmax','Tavg','Isp max','Isp avg','Impulse','t90%'})
title('Monte Carlo output spread')
ylabel('Value')

% Sensitivity-style scatter plots

if exist('aVec','var') && exist('OFVec','var') && exist('AtVec','var')
    figure('Name','Monte Carlo - Input/Output Sensitivity');

    subplot(2,3,1)
    scatter(OFVec, T_max, 12, 'filled')
    xlabel('O/F')
    ylabel('T_{max}')
    title('T_{max} vs O/F')
    grid on

    subplot(2,3,2)
    scatter(AtVec, T_max, 12, 'filled')
    xlabel('A_t')
    ylabel('T_{max}')
    title('T_{max} vs A_t')
    grid on

    subplot(2,3,3)
    scatter(AeVec, Isp_avg, 12, 'filled')
    xlabel('A_e')
    ylabel('Isp_{avg}')
    title('Isp_{avg} vs A_e')
    grid on

    subplot(2,3,4)
    scatter(pambVec, T_avg, 12, 'filled')
    xlabel('p_{amb}')
    ylabel('T_{avg}')
    title('T_{avg} vs p_{amb}')
    grid on

    subplot(2,3,5)
    scatter(dintVec, deltatGuaranteedThrust, 12, 'filled')
    xlabel('d_{int,0}')
    ylabel('t_{90\%}')
    title('t_{90\%} vs d_{int,0}')
    grid on

    subplot(2,3,6)
    scatter(aVec, rB, 12, 'filled')
    xlabel('a')
    ylabel('r_B')
    title('r_B vs a')
    grid on

    sgtitle('Monte Carlo Sensitivity Scatter Plots')
end

% Correlation matrix of main outputs
Y = [T_max(:), T_avg(:), Isp_max(:), Isp_avg(:), Isp_tot_max(:), deltatGuaranteedThrust(:)];
R = corrcoef(Y, 'Rows', 'complete');

figure('Name','Monte Carlo - Output Correlation');
imagesc(R)
colorbar
axis equal tight
set(gca, 'XTick', 1:6, 'XTickLabel', {'Tmax','Tavg','Isp max','Isp avg','Impulse','t90%'})
set(gca, 'YTick', 1:6, 'YTickLabel', {'Tmax','Tavg','Isp max','Isp avg','Impulse','t90%'})
title('Correlation matrix of Monte Carlo outputs')

% Normal probability plots (optional but useful)
figure('Name','Monte Carlo - Normality Check');

subplot(2,2,1)
qqplot(T_max)
title('QQ plot - T_{max}')

subplot(2,2,2)
qqplot(T_avg)
title('QQ plot - T_{avg}')

subplot(2,2,3)
qqplot(Isp_max)
title('QQ plot - Isp_{max}')

subplot(2,2,4)
qqplot(deltatGuaranteedThrust)
title('QQ plot - t_{90%}')