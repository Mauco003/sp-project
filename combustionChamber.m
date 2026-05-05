function results = combustionChamber(casing, liner, nozzle, propellant, performance, constants, options)
%% -----------------------------------------------------------------------
% Function to size the combustion chamber and characterize the heat
% transfer throughout it. 
%
% INPUTS:
% liner - struct including the liner parameters
%       .regressionRate
%       .thickness
%       .thermalConductivity
%
% OUTPUTS:
% results of thicknesses
% UNITS:
%
%
% EXAMPLE USAGE:
%-------------------------------------------------------------------------
%% Definition of variables
arguments
    casing                    struct
    liner                     struct
    nozzle                    struct
    propellant                struct
    performance               struct
    constants                 Constants
    options.makePlot    (1,1) logical = true
    options.showSummary (1,1) logical = true
end

% initial parameters for heat transfer
T_cc = propellant.cea.ccTemperature;     % [K] temperature of combustion chamber
%t_liner = liner.thickness;               % [mm] liner initial thickness
k_liner = liner.thermalConductivity;

r_burn = 8.88;                           % [mm/s]
t_burn = 25;                             % [s]

t_casing = casing.thickness;             % [mm]
k_casing = casing.thermalConductivity;
TMax = casing.TMax;

h_out = 10;                              % [W/m^2/K], natural convection
T_ambient = 293.15;                      % [K], assumed

% chemistry and gas composition
R = constants.R / propellant.cea.molarMass;

gamma = propellant.cea.gamma;
%pc    = propellant.ccPressure;
pc = 70e5;  
cStar = performance.cstar;
muGas = propellant.cea.mu;
cpGas = R * gamma/(gamma-1);
PrGas = muGas*cpGas / propellant.cea.k; % Suppose it is constant for now

% geometry of grain
d_grain = 0.89;
r_grain = d_grain/2;
l_grain = 1.58;

% geometry of nozzle
At    = nozzle.At;
alpha = nozzle.alpha;     % [rad]
beta  = nozzle.beta;      % [rad]
rt = sqrt(At/pi);
Dt = rt*2;

% use bartz correleation to estimate h
M = 0;  % in cc
T0 = T_cc;
A = pi*r_grain^2;
epsilon = A/At;
omega = 0.6;    % viscosity exponent in Bartz correction [-]
rCurve = nozzle.rCurvature;


% iterate to find the real Tw and hg
tol = 0.1;              % [K]
error = 100000;
error2 = 100000;
Tw_guess = 800;         %[K]
maxIter = 10000;
iter = 0;
liner_guess = 0.004;

% Outer-side total resistance
%R_out = t_liner/k_liner + t_casing/k_casing + 1/h_out; 

while((error > tol) && (iter < maxIter) && (error2 > tol))
    % using MEOP as pc.
    hg = bartzCorrelation(pc, cStar, Dt, rCurve, epsilon, muGas, cpGas, PrGas, Tw_guess, T_cc, T_cc, omega);

    R_out = liner_guess/k_liner + 1/hg;

    q = R_out*(T_cc - TMax);

    Tw_upd = T_cc - q/hg;

    t_liner = k_liner/q * abs(TMax - Tw_upd);

    % check error
    error = abs(Tw_upd - Tw_guess);
    error2 = abs(t_liner - liner_guess);

    % updating Tw guess dependant on the error
    
    Tw_guess = Tw_upd;
    liner_guess = t_liner;
    iter = iter + 1;
end

% % Heat flux
% q = hg * (T_cc - Tw_upd);   % [W/m^2]
% 
% %% Interface temperatures
% % 1 - convection from the gas to the liner
% T_linerInner = Tw_guess;
% 
% % 2 - conduction from the liner to the casing
% T_linerCasing = T_linerInner - q * (t_liner/k_liner);
% 
% % 3 - casing near side to far side ('to external environment')
% T_casingOuter = T_linerCasing - q * (t_casing/k_casing);


%% Sizing and optimization

% First check casing temperature
% if T_linerCasing < TMax
%     casing.survival = true;
% else
%     % we no gucci
%     casing.survival = false;
% end

% Next check lining exists for the entire burn time
linerConsumed = (liner.regressionRate* t_burn);   % [m]
disp(linerConsumed);
liner.survival = t_liner > linerConsumed;

if liner.survival
    % liner survives
    liner.survival = true;
else
    % liner is burned through before end of burn
    liner.survival = false;
end

% Results
results = struct();

% results.T_linerInner = T_linerInner;
% results.T_linerCasing = T_linerCasing;
% results.T_casingOuter = T_casingOuter;
results.q = q;
results.hg = hg;

%results.casingSurvives = casing.survival;
results.linerSurvives = liner.survival;
results.lt = liner_guess;

end