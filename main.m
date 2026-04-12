clear; close all; clc

addpath(genpath("."))
%% Data - Setup


% Pressure [bar]
pressureData = [10.1; 9.7; 10.3; 30.2; 31.0; 29.8; 50.2; 51.0; ...
            50.3; 70.2; 69.0; 69.8; 91.2; 89.1; 89.0];

% Burning rate [mm/s]
burningRateData = [4.0; 3.8; 4.1; 5.6; 6.0; 5.7; 7.0; 7.2; ...
                7.1; 8.4; 8.3; 8.6; 8.8; 9.0; 9.2];

[a, aSigma, n, nSigma, R2] = Uncertainty(pressureData, burningRateData);



%%
term1 = (2 * gamma^2) / (gamma - 1);
term2 = (2 / (gamma + 1))^((gamma + 1) / (gamma - 1));
term3 = 1 - (Pa / Pc)^((gamma - 1) / gamma);

Ct = sqrt(term1 * term2 * term3); % Thrust Coeff.
Pc_pa = Pc * 1e5; % Need Pascal for At    
At = T / (Pc_pa * Ct);

