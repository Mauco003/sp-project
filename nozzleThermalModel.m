function best = nozzleThermalModel(thermo, nozzle, cooling, options)
%NOZZLETHERMALMODEL
% Thermal model with:
%   - gas-side Bartz correlation
%   - one single conductive wall
%   - equivalent annular cooling jacket
%   - iterative search for minimum wall thickness and minimum water mass flow
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
    thermo struct
    nozzle struct
    cooling struct
    options.Npts (1,1) double = 10
    options.AA_limit (1,1) double = 2
    options.t_start (1,1) double = 1e-3      % 1 mm
    options.t_step  (1,1) double = 0.25e-3   % 0.25 mm
    options.t_max   (1,1) double = 20e-3
    options.mdotTol (1,1) double = 1e-4
    options.makePlot (1,1) logical = true
    options.showSummary (1,1) logical = true
end

%% ------------------------------------------------------------------------
% Fixed water properties (first version)
% -------------------------------------------------------------------------
cpWater  = 4180;        % [J/kg/K]
rhoWater = 997;         % [kg/m^3]
muWater  = 0.89e-3;     % [Pa*s]
kWater   = 0.60;        % [W/m/K]
PrWater  = cpWater * muWater / kWater; % Prandtl number

% inputs from previous work. I put structs bc there is way too many

% chemical

gamma = thermo.gamma;
Tc    = thermo.Tc;
pc    = thermo.pc;
cStar = thermo.cstar;
muGas = thermo.mu_gas;
cpGas = thermo.cp_gas;
PrGas = thermo.Pr_gas;

% nozzle

At    = nozzle.At;
alpha = nozzle.alpha;     % [rad]
beta  = nozzle.beta;      % [rad]
burnTime = nozzle.burnTime;

% Cooling jacket related
TwaterInitial = cooling.Twater_in;
kNozzleWall   = cooling.k_wall;
pWater        = cooling.pWater;
uWaterTarget  = cooling.uWaterTarget;

% Mesh for IT transfer

Npts = options.Npts;
areaRatioCool = options.AA_limit;

% Geometry of cooling jacket
% FIXED from At, alpha, beta and the desired area ratio

rt = sqrt(At/pi);
rLim = sqrt(areaRatioCool*At/pi);

Lconv = (rLim - rt)/tan(beta);
Ldiv  = (rLim - rt)/tan(alpha);

x = linspace(-Lconv, Ldiv, Npts).';
r = zeros(size(x));

for i = 1:Npts
    if x(i) <= 0
        r(i) = rt - x(i)*tan(beta);
    else
        r(i) = rt + x(i)*tan(alpha);
    end
end

% cross sectional area and Area ratio of each point of the mesh. Used to
% evaluate isentropic relationships

A = pi*r.^2;
AA = A/At;

% Hot flow quantities along the mesh. Initialized

M        = zeros(Npts,1);
T_static = zeros(Npts,1);
Taw      = zeros(Npts,1); % Adiabatic wall temperature
hg       = zeros(Npts,1); % COnvective coefficient for the core flow CHANGES

% Get the Mach number based on area ratio at each point of the mesh

for i = 1:Npts
    if abs(AA(i)-1) < 1e-8
        M(i) = 1.0;
    elseif x(i) < 0
        M(i) = machFromAreaRatio(AA(i), gamma, 'subsonic');
    else
        M(i) = machFromAreaRatio(AA(i), gamma, 'supersonic');
    end
    
    % Get temperature based on isentropic relationships

    T_static(i) = Tc * T_T0_Isen(M(i), gamma);
    Taw(i) = adiabaticWallTemperature(T_static(i), M(i), gamma, PrGas, thermo.recoveryModel);

    % Bartz correlation used for h
    Tw_ref = 1500;   % only for sigma correction? Sutton stuff
    hg(i) = bartzCorrelation(pc, cStar, 2*rt, A(i), At, ...
                             muGas, cpGas, PrGas, gamma, M(i), Tc, Tw_ref);
end

% Boiling temperature from water pressure

TBoil = waterSaturationTemperature(pWater); % depends on pressure!

% Thickness of the wall loop to make it as thin as possible without water
% boiling

% FIX THE TWALL. THE SAME? FOR ALL NOZZLE? Or should it be according the T
% adiabatic

% consider a regession rate for TBC
% dont model it

%% ------------------------------------------------------------------------
% Summary
% -------------------------------------------------------------------------
if options.showSummary
    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('              NOZZLE THERMAL MODEL                  \n');
    fprintf('====================================================\n');
    fprintf('Selected wall thickness      : %.6f m\n', best.t_wall);
    fprintf('Minimum water mass flow      : %.6f kg/s\n', best.cooling.mDotWater);
    fprintf('Equivalent channel gap       : %.6f m\n', best.cooling.gap);
    fprintf('Equivalent flow area         : %.6e m^2\n', best.cooling.Aflow);
    fprintf('Hydraulic diameter           : %.6f m\n', best.cooling.Dh);
    fprintf('Water Reynolds               : %.3f\n', best.cooling.Re);
    fprintf('Water Nusselt                : %.3f\n', best.cooling.Nu);
    fprintf('Water-side h                 : %.3f W/m^2/K\n', best.cooling.hWater);
    fprintf('Coolant inlet temperature    : %.3f K\n', TwaterInitial);
    fprintf('Coolant outlet temperature   : %.3f K\n', best.Twater_out);
    fprintf('Boiling limit                : %.3f K\n', best.cooling.TBoil);
    fprintf('Total heat rate              : %.3f W\n', best.Qdot_total);
    fprintf('Total heat                   : %.3f J\n', best.Q_total);
    fprintf('====================================================\n');
    fprintf('\n');
end

%% ------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
if options.makePlot
    figure;

    subplot(2,2,1)
    plot(x, T_static, 'LineWidth', 1.5); hold on
    plot(x, Taw, 'LineWidth', 1.5);
    plot(x, best.Tw_hot, 'LineWidth', 1.5);
    plot(x, best.Tw_cold, 'LineWidth', 1.5);
    plot(x, best.Twater, 'LineWidth', 1.5);
    grid on
    xlabel('x [m]')
    ylabel('Temperature [K]')
    legend('T_{static}','T_{aw}','T_{w,hot}','T_{w,cold}','T_{water}','Location','best')
    title(sprintf('Thermal profile, t = %.2f mm', 1e3*best.t_wall))

    subplot(2,2,2)
    plot(x, M, 'LineWidth', 1.5)
    grid on
    xlabel('x [m]')
    ylabel('Mach')
    title('Mach profile')

    subplot(2,2,3)
    plot(x, hg, 'LineWidth', 1.5)
    grid on
    xlabel('x [m]')
    ylabel('h_g [W/m^2/K]')
    title('Bartz coefficient')

    subplot(2,2,4)
    xmid = 0.5*(x(1:end-1)+x(2:end));
    plot(xmid, best.qpp, 'LineWidth', 1.5)
    grid on
    xlabel('x [m]')
    ylabel('q'''' [W/m^2]')
    title('Heat flux')
end

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