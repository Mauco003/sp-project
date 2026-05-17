function results = combustionChamber(casing, liner, nozzle, grain, propellant, performance, constants)
%% -----------------------------------------------------------------------
% Function to size the combustion chamber and nozzle's ablative materials, and characterize the heat
% transfer throughout it. 
%
% INPUTS:
% casing - struct including the casing parameters
%       .thickness
%       .density
%       .hoopStress
%       .safetyFactor
%       .thermalConductivity
%       .TMax
%       .cost
% liner - struct including the liner parameters
%       .regressionRate
%       .thickness
%       .thermalConductivity
%       .density
%       .cost
% nozzle - struct including the nozzle parameters
%       .At
%       .Ae
%       .alpha
%       .beta
%       .rCurvature
% grain - struct including the grain parameters
%       .dInt0
%       .dExt0
%       .L0
%
% propellant - struct including the propellant parameters
%       .cea - struct including the CEA outputs
%           .ccTemperature
%           .gamma
%           .molarMass
%           .mu
%           .k
%           .ccPressure
% performance - struct including the performance parameters
%       .cstar
%       .mDot
%
% OUTPUTS:
% results - struct including the results of the thermal model and optimization
%       .q - heat flux to the wall [W/m^2]
%       .hg - convective coefficient at the wall [W/m^2/K]
%       .linerSurvives - boolean indicating whether the liner survives the burn time
%       .casingSurvives - boolean indicating whether the casing survives the burn time
%       .lt - liner thickness required for survival [m]
%
%
% EXAMPLE USAGE:
% results = combustionChamber(casing, liner, nozzle, propellant, performance, constants);
%-------------------------------------------------------------------------
%% Definition of variables
arguments
    casing                    struct
    liner                     struct
    nozzle                    struct
    grain                     struct
    propellant                struct
    performance               struct
    constants                 Constants
end

% Getting input data
T_cc = propellant.cea.ccTemperature;     % [K] temperature of combustion chamber
k_liner = liner.thermalConductivity;

% Constants
r_burn = constants.rb;                                % [mm/s]
t_burn = constants.tBurn;                             % [s]
h_out = constants.h_out;                              % [W/m^2/K], natural convection
T_ambient = constants.T_ambient;                      % [K], assumed

% Casing, liner, and other parameters
t_casing = casing.thickness;                          % [mm]
k_casing = casing.thermalConductivity;
TMax = casing.TMax;

% Grain geometry
d_grain = grain.dInt0;
r_grain = d_grain/2;
l_grain = grain.L0;

% Nozzle geometry
At    = nozzle.At;
alpha = nozzle.alpha;     % [rad]
beta  = nozzle.beta;      % [rad]
rt = sqrt(At/pi);
Dt = rt*2;
rCurve = nozzle.rCurvature;

% Chemistry and gas composition
R = constants.R / propellant.cea.molarMass;
gamma = propellant.cea.gamma;
pc    = propellant.ccPressure;
cStar = performance.cstar;
muGas = propellant.cea.mu;
cpGas = R * gamma/(gamma-1);
PrGas = muGas*cpGas / propellant.cea.k; 


% Use bartz correleation to estimate h
A = pi*r_grain^2;
epsilon = A/At;
omega = 0.6;    % viscosity exponent in Bartz correction [-]


% Itertion parameters to find the real Tw and hg
tol = 0.1;              % [K]
error = 100000;
error2 = 100000;
Tw_guess = 800;         %[K]
maxIter = 10000;
iter = 0;
liner_guess = 0.004;

% Iteration loop 
while((error > tol) && (iter < maxIter) && (error2 > tol))
    % using MEOP as pc.
    hg = bartzCorrelation(pc(1), cStar, Dt, rCurve, epsilon, muGas, cpGas, PrGas, Tw_guess, T_cc, T_cc, omega);

    % Compute heat transfer problem
    R_out = liner_guess/k_liner + 1/hg;
    q = R_out*(T_cc - TMax);
    Tw_upd = T_cc - q/hg;

    % Thickness of the liner based on Tw
    t_liner = k_liner/q * abs(TMax - Tw_upd);

    % Update errors
    error = abs(Tw_upd - Tw_guess);
    error2 = abs(t_liner - liner_guess);

    % Updating Tw guess dependant on the error
    Tw_guess = Tw_upd;
    liner_guess = t_liner;
    iter = iter + 1;

end

%% Sizing and optimization
% Next check liner exists for the entire burn time
linerConsumed = (liner.regressionRate* t_burn);   % [m]

% Total liner thickness required:
% thermal protection thickness + consumed ablative thickness
linerRequired = liner_guess + linerConsumed;     % [m]

% Check if selected candidate liner thickness is sufficient
liner.survival = liner.thickness >= linerRequired;

% Check casing hoop stress
sigmaHoop = pc(1) * grain.dExt0 / (2 * casing.thickness);   % [Pa]

% Allowable stress including safety factor
sigmaAllow = casing.hoopStress / casing.safetyFactor;        % [Pa]

casing.survival = sigmaHoop <= sigmaAllow;


%% Results
results = struct();

results.q = q;
results.hg = hg;
results.casingSurvives = casing.survival;
results.linerSurvives = liner.survival;
results.liner_thickness = liner_guess;
results.linerRequired = linerRequired;
results.thermalMargin = liner.thickness / linerRequired;
results.structuralMargin = sigmaAllow / sigmaHoop;

end