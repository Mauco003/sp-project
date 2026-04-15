function best = nozzleThermalModel(input, propellant, nozzle, performance, cooling, constants, options)
%NOZZLETHERMALMODEL
% Thermal model with:
%   - gas-side Bartz correlation
%   - one single conductive wall
%   - equivalent annular cooling jacket
%   - iterative search for minimum wall thickness and minimum water mass flow
%
% Conical nozzle is assumed
%
% Main idea:
%   q'' = (Taw - Twater_bulk) / (1/hg + t/k + 1/hWater)
%
% Inputs:
%   thermo.gamma
%   thermo.Tc
%   thermo.pc
%   thermo.cstar
%   thermo.mu_gas
%   thermo.cp_gas
%   thermo.Pr_gas
%   thermo.recoveryModel
%
%   nozzle.At
%   nozzle.alpha
%   nozzle.beta
%   nozzle.burnTime
%
%   cooling.Twater_in
%   cooling.k_wall
%   cooling.pWater         % coolant pressure [Pa]
%   cooling.uWaterTarget   % target bulk velocity in cooling jacket [m/s]
%   cooling.TBCthickness   % [m]
%
% Options:
%   options.Npts
%   options.AA_limit
%   options.t_start
%   options.t_step
%   options.t_max
%   options.mdotTol
%   options.makePlot
%   options.showSummary

arguments
    input
    propellant                struct
    nozzle                    struct
    performance               struct
    cooling                   struct
    constants                 Constants
    options.n           (1,1) double  = 10        % [-] Number of points for each section (converging section), total points = 2n - 1
    options.areaLimits  (2,1) double  = [2; 2]    % The jacket goes from A(1) to A(2)
    options.t_start     (1,1) double  = 1e-3      % 1 mm
    options.t_step      (1,1) double  = 0.25e-3   % 0.25 mm
    options.t_max       (1,1) double  = 20e-3
    options.mdotTol     (1,1) double  = 1e-4
    options.makePlot    (1,1) logical = true
    options.showSummary (1,1) logical = true
end

%% ------------------------------------------------------------------------
% Fixed water properties (first version)
% -------------------------------------------------------------------------
% cooling.cp  = 4180;                                 % [J/kg/K]
% cooling.rho = 997;                                  % [kg/m^3]
% cooling.mu  = 0.89e-3;                              % [Pa*s]
% cooling.k   = 0.60;                                 % [W/m/K]
% cooling.Pr  = cooling.cp * cooling.mu / cooling.k;  % Prandtl number

waterT0         = cooling(end).temperature;
% kNozzleWall   = cooling.kWall;
waterPressure   = cooling.pressure;
% uWaterTarget  = cooling.velocity;

% inputs from previous work. I put structs bc there is way too many

% chemical
R = constants.R / propellant.cea.molarMass;

gamma = propellant.cea.gamma;
Tc    = propellant.cea.ccTemperature;
pc    = propellant.ccPressure;
cStar = performance.cstar;
muGas = propellant.cea.mu;
cpGas = R * gamma/(gamma-1);
PrGas = muGas*cpGas / propellant.cea.k;

% nozzle
At    = nozzle.At;
alpha = nozzle.alpha;     % [rad]
beta  = nozzle.beta;      % [rad]
% burnTime = performance.t; % [s]

% Mesh for IT transfer
n = options.n;
areaRatioCool = options.areaLimits;

% Geometry of cooling jacket
% FIXED from At, alpha, beta and the desired area ratio
rt = sqrt(At/pi);
rLim = sqrt(areaRatioCool*At/pi);

Lconv = (rLim(1) - rt)/tan(beta);
Ldiv  = (rLim(2) - rt)/tan(alpha);

xConv = linspace(-Lconv, 0, n)';
xDiv  = linspace(0, Ldiv, n)';

x = [xConv; xDiv(2:end)];
r = zeros(size(x));

r(x <= 0) = rt - x(x <= 0)*tan(beta);
r(x > 0) = rt + x(x > 0)*tan(alpha);

% cross sectional area and Area ratio of each point of the mesh. Used to
% evaluate isentropic relationships
A = pi*r.^2;
epsilon = A/At;

% Hot flow quantities along the mesh. Initialized
% Get the Mach number based on area ratio at each point of the mesh
mach        = zeros(2*n-1,1);
mach(n) = 1;
for i = 1:n-1,         mach(i) = machFromAreaRatio(epsilon(i), gamma, 'subsonic');   end
for i = (n+1):(2*n-1), mach(i) = machFromAreaRatio(epsilon(i), gamma, 'supersonic'); end

T0 = Tc * (1 + 0.5*(gamma-1)*performance.Mcc^2); % Under adiabatic conditions, T0 is constant

% Bartz correlation used for h
Twall = 1500;   % only for sigma correction? Sutton stuff

% Computing heat flux
% Data i need:
%   T(x), isoentropic
%   T, adiabatic

% Table of temperature along mesh for each cooling point
T = zeros(2*n-1, length(cooling) + 1);
dx = zeros(2*n-1,2); % Thickness of each wall section


for i = 1:length(x)
    % Solve heat flux at station i
    recoveryFactor = (1 + 0.5*PrGas^(1/3)*(gamma-1)*mach.^2)./ ...
    (1+ 0.5*(gamma-1)*mach.^2);                  % Based on 05-Liquid-PART4-Cooling, pg. 10
    Taw = T0 .*recoveryFactor;
    
    if i == 1
        T(i, end) = cooling(end).temperature; % Initial guess for water temperature at inlet
        T(i, 1) = Taw;
    else
        T(i, 1) = T(i-1, 1);
        T(i, end) = T(i-1, end) + computeHeatFlux(cooling, x(i-1:i), r(i-1:i), T(i-1:i, :));
    end


    h1 = bartzCorrelation(pc, cstar, 2*rt, nozzle.rCurvature, epsilon, mu, cp, Pr);
    k2 = cooling(2).k;
    k3 = cooling(3).k;
    h4 = 1; % Celia

    H = 1/(1/h1 + dx(i, 1)/k2 + dx(i, 2)/k3 + 1/h4);

    q = H * (T(i, 1) - T(i, end));




end

[q] = computeHeatFlux(cooling, x, r);

% TBoil = waterSaturationTemperature(pWater); % depends on pressure!

% Thickness of the wall loop to make it as thin as possible without water
% boiling

% FIX THE TWALL. THE SAME? FOR ALL NOZZLE? Or should it be according the T
% adiabatic

% consider a regession rate for TBC
% dont model it

%% ------------------------------------------------------------------------
% Summary
% -------------------------------------------------------------------------
% if options.showSummary
%     fprintf('\n');
%     fprintf('====================================================\n');
%     fprintf('              NOZZLE THERMAL MODEL                  \n');
%     fprintf('====================================================\n');
%     fprintf('Selected wall thickness      : %.6f m\n', best.t_wall);
%     fprintf('Minimum water mass flow      : %.6f kg/s\n', best.cooling.mDotWater);
%     fprintf('Equivalent channel gap       : %.6f m\n', best.cooling.gap);
%     fprintf('Equivalent flow area         : %.6e m^2\n', best.cooling.Aflow);
%     fprintf('Hydraulic diameter           : %.6f m\n', best.cooling.Dh);
%     fprintf('Water Reynolds               : %.3f\n', best.cooling.Re);
%     fprintf('Water Nusselt                : %.3f\n', best.cooling.Nu);
%     fprintf('Water-side h                 : %.3f W/m^2/K\n', best.cooling.hWater);
%     fprintf('Coolant inlet temperature    : %.3f K\n', TwaterInitial);
%     fprintf('Coolant outlet temperature   : %.3f K\n', best.Twater_out);
%     fprintf('Boiling limit                : %.3f K\n', best.cooling.TBoil);
%     fprintf('Total heat rate              : %.3f W\n', best.Qdot_total);
%     fprintf('Total heat                   : %.3f J\n', best.Q_total);
%     fprintf('====================================================\n');
%     fprintf('\n');
% end

%% ------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
% if options.makePlot
%     figure;
% 
%     subplot(2,2,1)
%     plot(x, T0, 'LineWidth', 1.5); hold on
%     plot(x, Taw, 'LineWidth', 1.5);
%     plot(x, best.Tw_hot, 'LineWidth', 1.5);
%     plot(x, best.Tw_cold, 'LineWidth', 1.5);
%     plot(x, best.Twater, 'LineWidth', 1.5);
%     grid on
%     xlabel('x [m]')
%     ylabel('Temperature [K]')
%     legend('T_{static}','T_{aw}','T_{w,hot}','T_{w,cold}','T_{water}','Location','best')
%     title(sprintf('Thermal profile, t = %.2f mm', 1e3*best.t_wall))
% 
%     subplot(2,2,2)
%     plot(x, mach, 'LineWidth', 1.5)
%     grid on
%     xlabel('x [m]')
%     ylabel('Mach')
%     title('Mach profile')
% 
%     subplot(2,2,3)
%     plot(x, hg, 'LineWidth', 1.5)
%     grid on
%     xlabel('x [m]')
%     ylabel('h_g [W/m^2/K]')
%     title('Bartz coefficient')
% 
%     subplot(2,2,4)
%     xmid = 0.5*(x(1:end-1)+x(2:end));
%     plot(xmid, best.qpp, 'LineWidth', 1.5)
%     grid on
%     xlabel('x [m]')
%     ylabel('q'''' [W/m^2]')
%     title('Heat flux')
% end

end


%% ========================================================================
% Helper: simulation for one wall thickness and one water mass flow
% ========================================================================
function result = simulateOneCase(mDotWater, tWall, x, r, A, AA, ...
                                  M, T_static, Taw, hg, ...
                                  rhoWater, muWater, kWater, PrWater, cpWater, ...
                                  TwaterInitial, uWaterTarget, kNozzleWall)

Npts = length(x);

Pmean = mean(2*pi*r);
Aflow = mDotWater / (rhoWater * uWaterTarget);
gap   = Aflow / Pmean;   % equivalent annular gap
Dh    = 2*gap;           % thin annulus approximation

Re = rhoWater * uWaterTarget * Dh / muWater;

if Re < 2300
    Nu = 4.36;   % simple laminar fully developed approximation
else
    Nu = 0.023 * Re^0.8 * PrWater^0.4;  % Dittus-Boelter
end

hWater = Nu * kWater / Dh;

Twater  = zeros(Npts,1);
Tw_hot  = zeros(Npts,1);
Tw_cold = zeros(Npts,1);

qpp      = zeros(Npts-1,1);
Qdot_seg = zeros(Npts-1,1);
dA_seg   = zeros(Npts-1,1);

Twater(1) = TwaterInitial;

for i = 1:Npts

    Rtot = 1/hg(i) + tWall/kNozzleWall + 1/hWater;

    qpp_i = (Taw(i) - Twater(i)) / Rtot;

    Tw_hot(i)  = Taw(i) - qpp_i/hg(i);
    Tw_cold(i) = Tw_hot(i) - qpp_i*(tWall/kNozzleWall);

    if i < Npts
        dA_seg(i) = frustumStripArea(x(i), r(i), x(i+1), r(i+1));
        qpp(i)    = qpp_i;
        Qdot_seg(i) = qpp(i) * dA_seg(i);

        Twater(i+1) = Twater(i) + Qdot_seg(i)/(mDotWater*cpWater);
    end
end

result.Twater = Twater;
result.Twater_out = Twater(end);

result.Tw_hot = Tw_hot;
result.Tw_cold = Tw_cold;

result.qpp = qpp;
result.Qdot_seg = Qdot_seg;

result.Aflow = Aflow;
result.gap = gap;
result.Dh = Dh;
result.Re = Re;
result.Nu = Nu;
result.hWater = hWater;
end


%% ========================================================================
% Helper: saturation temperature of water from pressure
% Antoine correlation, rough engineering use
% ========================================================================
function Tsat = waterSaturationTemperature(pPa)
% pressure in Pa -> Tsat in K

p_mmHg = pPa / 133.322;

% Antoine constants for water (rough engineering range)
A = 8.14019;
B = 1810.94;
C = 244.485;

T_C = B / (A - log10(p_mmHg)) - C;
Tsat = T_C + 273.15;
end