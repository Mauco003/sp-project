function qDotAlongX = nozzleThermalModel(input, propellant, nozzle, performance, cooling, constants, options)

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
    options.compareCEA   (1,1) logical = true
end

cooling.cp  = 4180;                                 % [J/kg/K]
cooling.rho = 997;                                  % [kg/m^3]
cooling.mu  = 0.89e-3;                              % [Pa*s]
cooling.k   = 0.60;                                 % [W/m/K]
cooling.Pr  = cooling.cp * cooling.mu / cooling.k;  % Prandtl number

waterT0       = cooling(end).temperature;
kNozzleWall   = cooling.kWall;
waterPressure = cooling.pressure;

R = constants.R / propellant.cea.molarMass;

gamma = propellant.cea.gamma;
Tc    = propellant.cea.ccTemperature;
pc    = propellant.ccPressure;
cStar = performance.cstar;
muGas = propellant.cea.mu;
cpGas = R * gamma/(gamma-1);
PrGas = muGas*cpGas / propellant.cea.k; % Suppose it is constant for now

% nozzle
At    = nozzle.At;
alpha = nozzle.alpha;     % [rad]
beta  = nozzle.beta;      % [rad]
% burnTime = performance.t; % [s]

% Mesh for IT transfer
n = options.n;
N = 2*n - 1; % Total number of points in mesh
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

% Get data from CEA

if options.compareCEA
    epsConvCEA = sort(epsilon(1:n), 'descend')';
    epsDivCEA  = sort(epsilon(n+1:end), 'ascend')';

    % adjust BC if needed
    [~, GasCEA] = getThermoProfileCEA(80, 20, input.pcNominal/1e5, epsConvCEA, epsDivCEA);

    % Convergent branch of nozzle mesh
    epsConvMesh = epsilon(1:n);

    % Divergent branch of nozzle mesh
    epsDivMesh  = epsilon(n:N);

    % CEA data branches
    epsConvCEA_data   = GasCEA.eps(1:length(epsConvCEA));
    MachConvCEA_data  = GasCEA.Mach(1:length(epsConvCEA));
    TConvCEA_data     = GasCEA.T(1:length(epsConvCEA));
    gammaConvCEA_data = GasCEA.gamma(1:length(epsConvCEA));
    cpConvCEA_data    = GasCEA.cp(1:length(epsConvCEA));

    % throat + divergent branch
    epsDivCEA_data   = [1, GasCEA.eps(length(epsConvCEA)+2:end)];
    MachDivCEA_data  = [GasCEA.Mach(length(epsConvCEA)+1),  GasCEA.Mach(length(epsConvCEA)+2:end)];
    TDivCEA_data     = [GasCEA.T(length(epsConvCEA)+1),     GasCEA.T(length(epsConvCEA)+2:end)];
    gammaDivCEA_data = [GasCEA.gamma(length(epsConvCEA)+1), GasCEA.gamma(length(epsConvCEA)+2:end)];
    cpDivCEA_data    = [GasCEA.cp(length(epsConvCEA)+1),    GasCEA.cp(length(epsConvCEA)+2:end)];

    % Interpolate each branch separately because i have problems otherwise
    machCEA  = zeros(size(epsilon));
    TstatCEA = zeros(size(epsilon));
    gammaCEA = zeros(size(epsilon));
    cpCEA    = zeros(size(epsilon));

    machCEA(1:n)  = interp1(epsConvCEA_data, MachConvCEA_data, epsConvMesh, 'linear', 'extrap');
    TstatCEA(1:n) = interp1(epsConvCEA_data, TConvCEA_data,    epsConvMesh, 'linear', 'extrap');
    gammaCEA(1:n) = interp1(epsConvCEA_data, gammaConvCEA_data, epsConvMesh, 'linear', 'extrap');
    cpCEA(1:n)    = interp1(epsConvCEA_data, cpConvCEA_data,   epsConvMesh, 'linear', 'extrap');

    machCEA(n:N)  = interp1(epsDivCEA_data, MachDivCEA_data, epsDivMesh, 'linear', 'extrap');
    TstatCEA(n:N) = interp1(epsDivCEA_data, TDivCEA_data,    epsDivMesh, 'linear', 'extrap');
    gammaCEA(n:N) = interp1(epsDivCEA_data, gammaDivCEA_data, epsDivMesh, 'linear', 'extrap');
    cpCEA(n:N)    = interp1(epsDivCEA_data, cpDivCEA_data,   epsDivMesh, 'linear', 'extrap');

    % First approximation: keep mu and Pr constant
    muCEA = muGas * ones(size(epsilon));
    PrCEA = PrGas * ones(size(epsilon));
else
    machCEA  = [];
    TstatCEA = [];
    gammaCEA = [];
    cpCEA    = [];
    muCEA    = [];
    PrCEA    = [];
end

%% ISENTROPIC APPROACH

% Get the Mach number based on area ratio at each point of the mesh
mach        = zeros(N,1);
mach(n) = 1;
for i = 1:n-1,         mach(i) = machFromAreaRatio(epsilon(i), gamma, 'subsonic');   end
for i = (n+1):(N), mach(i) = machFromAreaRatio(epsilon(i), gamma, 'supersonic'); end

T0 = Tc * (1 + 0.5*(gamma-1)*performance.Mcc^2); % Under adiabatic conditions, T0 is constant

% Bartz correlation used for h
% Twall = 1500;

% Computing heat flux
% Data i need:
%   T(x), isoentropic
%   T, adiabatic

% Table of temperature along mesh for each cooling point
T = zeros(N, length(cooling) + 1);

dx = zeros(N,2);
SF = 1.5; % Same as everything
dx(:,1) = cfg.tWall;   % wall thickness [m]
dx(:,2) = best.tTBC*SF;   % TBC thickness [m]
mDot = 3; % [kg/s] obtained for mas q_dot and a temp raise of 80 degrees in bulk as first approach
q  = zeros(N-1,1); % Heat flux at each section

T(1, end) = cooling(end).temperature; % Initial guess for water temperature at inlet

% Surface area of each segment of the mesh, used for heat flux calculation

dA = 2*pi*(sqrt(diff(x).^2 + diff(r).^2)) .* 0.5 .* (r(1:end-1) + r(2:end));

for i = 1:N
    T(i, 1) = T0*recoveryFactor(gamma, mach(i), PrGas); % Adiabatic wall temperature at station i

    h1 = bartzCorrelation(input.pcNominal, performance.cstar, 2*rt, nozzle.rCurvature, epsilon(i), muGas, cpGas, PrGas);
    k2 = cooling(2).k;
    k3 = cooling(3).k;
    h4 = dittusCorrelation(cooling(end).Dc,cooling(end).k,cooling(end).Re, cooling(end).Pr); %dittusCorrelation(Dc,kH2O,Re,Pr);

    H = 1/(1/h1 + dx(i, 1)/k2 + dx(i, 2)/k3 + 1/h4);

    q(i) = H * (T(i, 1) - T(i, end));

    % Update intermediate station tempeartures
    T(i, 2) = T(i, 1) - q(i)/h1;
    T(i, 3) = T(i, 2) - q(i)*dx(i, 1)/k2;
    T(i, 4) = T(i, 3) - q(i)*dx(i, 2)/k3;

    if i < N
        T(i+1, end) = T(i, end) + q(i)*dA(i)/(cooling(end).cp*cooling(end).mDot); % Update water temperature at station i
    end
end

%% CEA
if options.compareCEA

    T_CEA  = zeros(N, length(cooling) + 1);
    q_CEA  = zeros(N,1);
    h1_CEA = zeros(N,1);
    T0_CEA = zeros(N,1);

    T_CEA(1,end) = cooling(end).temperature;   % coolant inlet
    
    for i = 1:N
        % Adiabatic wall temperature using local CEA data
        T0_CEA(i) = TstatCEA(i) * (1 + 0.5*(gammaCEA(i)-1)*machCEA(i)^2); %%% THIS MIGHT BE WRONG. OBTAIN TOTAL TEMP FROM CEAAA
        T_CEA(i,1) = T0_CEA(i) * recoveryFactor(gammaCEA(i), machCEA(i), PrCEA(i));

        h1_CEA(i) = bartzCorrelation(input.pcNominal, performance.cstar, ...
            2*rt, nozzle.rCurvature, epsilon(i), muCEA(i), cpCEA(i), PrCEA(i));

        k2 = cooling(2).k;
        k3 = cooling(3).k;
        h4 = dittusCorrelation(cooling(end).Dc, cooling(end).k, ...
            cooling(end).Re, cooling(end).Pr);

        H_CEA = 1/(1/h1_CEA(i) + dx(i,1)/k2 + dx(i,2)/k3 + 1/h4);

        q_CEA(i) = H_CEA * (T_CEA(i,1) - T_CEA(i,end));

        T_CEA(i,2) = T_CEA(i,1) - q_CEA(i)/h1_CEA(i);
        T_CEA(i,3) = T_CEA(i,2) - q_CEA(i)*dx(i,1)/k2;
        T_CEA(i,4) = T_CEA(i,3) - q_CEA(i)*dx(i,2)/k3;

        if i < N
            T_CEA(i+1,end) = T_CEA(i,end) + q_CEA(i)*dA(i)/(cooling(end).cp*cooling(end).mDot);
        end
    end

else

    T_CEA  = [];
    q_CEA  = [];
    h1_CEA = [];
    T0_CEA = [];

end

best = struct();

best.x       = x;
best.r       = r;
best.A       = A;
best.epsilon = epsilon;
best.mach    = mach;

best.T       = T;
best.q       = q;
best.dA      = dA;
best.T0_CEA = T0_CEA;

best.dx      = dx;
best.waterInletTemperature = cooling(end).temperature;
best.waterPressure         = cooling.pressure;

best.T_CEA   = T_CEA;
best.q_CEA   = q_CEA;
best.h1_CEA  = h1_CEA;

best.machCEA  = machCEA;
best.TstatCEA = TstatCEA;
best.gammaCEA = gammaCEA;
best.cpCEA    = cpCEA;

%% PLOTS

if options.makePlot

    % Coolant temperature
    figure('Name','Coolant temperature','NumberTitle','off');
    plot(x, T(:,end) - 273.15, 'LineWidth', 1.8); hold on
    if options.compareCEA
        plot(x, T_CEA(:,end) - 273.15, '--', 'LineWidth', 1.8)
        legend('Ideal/perfect gas','CEA','Location','best')
    end
    grid on
    xlabel('x [m]')
    ylabel('Coolant temperature [^\circ C]')
    title('Coolant temperature along cooling jacket')

    % Wall temperature (hot side)
    figure('Name','Hot wall temperature','NumberTitle','off');
    plot(x, T(:,2), 'LineWidth', 1.8); hold on
    if options.compareCEA
        plot(x, T_CEA(:,2), '--', 'LineWidth', 1.8)
        legend('Ideal/perfect gas','CEA','Location','best')
    end
    grid on
    xlabel('x [m]')
    ylabel('Wall temperature [K]')
    title('Hot-side wall temperature')

    % Temperature through the wall / layers
    figure('Name','Temperature through layers','NumberTitle','off');
    plot(x, T(:,1), 'LineWidth', 1.8); hold on
    plot(x, T(:,2), 'LineWidth', 1.8)
    plot(x, T(:,3), 'LineWidth', 1.8)
    plot(x, T(:,4), 'LineWidth', 1.8)
    plot(x, T(:,end), 'LineWidth', 1.8)

    if options.compareCEA
        plot(x, T_CEA(:,1), '--', 'LineWidth', 1.8)
        plot(x, T_CEA(:,2), '--', 'LineWidth', 1.8)
        plot(x, T_CEA(:,3), '--', 'LineWidth', 1.8)
        plot(x, T_CEA(:,4), '--', 'LineWidth', 1.8)
        plot(x, T_CEA(:,end), '--', 'LineWidth', 1.8)

        legend('T_{aw} ideal','Wall hot ideal','After wall ideal','After TBC ideal','Coolant ideal', ...
               'T_{aw} CEA','Wall hot CEA','After wall CEA','After TBC CEA','Coolant CEA', ...
               'Location','best')
    else
        legend('T_{aw}','Wall hot side','After wall','After TBC','Coolant', ...
            'Location','best')
    end
    grid on
    xlabel('x [m]')
    ylabel('Temperature [K]')
    title('Temperature profile across layers')

    % Heat flux
    figure('Name','Heat flux','NumberTitle','off');
    plot(x, q, 'LineWidth', 1.8); hold on
    if options.compareCEA
        plot(x, q_CEA, '--', 'LineWidth', 1.8)
        legend('Ideal/perfect gas','CEA','Location','best')
    end
    grid on
    xlabel('x [m]')
    ylabel('Heat flux [W/m^2]')
    title('Heat flux along nozzle')

    figure;
    plot(x, machCEA, 'LineWidth', 1.8)
    grid on
    xlabel('x [m]')
    ylabel('Mach [-]')
    title('Mach from CEA interpolated to nozzle mesh')
    figure;
    plot(x, T0_CEA, 'LineWidth', 1.8)
    grid on
    xlabel('x [m]')
    ylabel('T_{CEA} [K]')
    title('Static temperature from CEA')

    figure;
    plot(x, gammaCEA, 'LineWidth', 1.8)
    grid on
    xlabel('x [m]')
    ylabel('\gamma_{CEA} [-]')
    title('Gamma from CEA')

    figure;
    plot(x, cpCEA, 'LineWidth', 1.8)
    grid on
    xlabel('x [m]')
    ylabel('c_p [J/kg/K]')
    title('cp from CEA')
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