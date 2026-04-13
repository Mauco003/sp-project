function nozzle = nozzleDesign(thrust, tc, pc, pe, gamma, molarMass, constants, options)
%NOZZLEDESIGN Preliminary nozzle geometry model
%
arguments
    thrust
    tc
    pc
    pe
    gamma
    molarMass
    constants   Constants = Constants()
    options.alpha = 15 * pi/180
    options.beta  = 30 * pi/180
    options.lambdaDiv = []
    options.machCC = 0.2
    options.ccRadius = []
    options.plot = false
    options.showSummary = false
end

%% Defaults / optional inputs
if ~isempty(options.lambdaDiv)
    alpha = acos(2*options.lambdaDiv - 1);
else
    alpha = options.alpha;
    lambdaDiv = (1 + cos(alpha))/2;
end

beta = options.beta;
machCC = options.machCC;

g0    = constants.g0;
R     = constants.R / molarMass;

% Ideal thermodynamic performance
ve    = exhaustVelocityIdeal(gamma, R, tc, pe, pc);
cstar = cstarIdeal(R, tc, gamma);
Isp   = ve / g0;
mdot  = thrust/(Isp*g0);

% Throat and exit areas
cf = cfIdeal(gamma, pe, pc);

At = thrust / (pc * cf);
epsilon = computeEpsilon(gamma, pe, pc);
Ae = epsilon * At;

rt = sqrt(At/pi);
re = sqrt(Ae/pi);

% Conical divergence length
lDiv = (re - rt) / tan(alpha);

% Chamber section from machCC if r_cc is not given
if ~isempty(options.ccRadius)
    rcc = ccRadius;
    Acc = pi * rcc^2;
    machCC = machFromAreaRatio(Acc/At, gamma, 'subsonic');   % not used because r_cc was directly imposed
else
    Acc = At * areaRatioIsen(machCC, gamma);
    rcc = sqrt(Acc/pi);
end

% Convergent section
lConv = (rcc - rt) / tan(beta);

% Geometry by nozzle type
lTotal = lConv + lDiv;

% Simple piecewise geometry for plotting
xConv = linspace(-lConv, 0, 100);
rConv = linspace(rcc, rt, 100);

xDiv = linspace(0, lDiv, 150);
rDiv = linspace(rt, re, 150);

x = [xConv, xDiv];
r = [rConv, rDiv];

% Exporting data
nozzle.x = [0, lConv, lTotal];
nozzle.r = [rcc, rt, re];

nozzle.Acc = Acc;
nozzle.At = At;
nozzle.Ae = Ae;

nozzle.ve    = ve;
% nozzle.cstar = cstar;
nozzle.Isp   = Isp;
nozzle.mDot  = mdot;
nozzle.cf    = cf;
nozzle.cstar = cstar;
nozzle.machCC = machCC;

% Plot
if options.plot
    figure;
    plot(x,  r, 'LineWidth', 1.5); hold on;
    plot(x, -r, 'LineWidth', 1.5);
    axis equal;
    grid on;
    xlabel('x [m]');
    ylabel('r [m]');
    title('Axisymmetric conical nozzle');
end

% Print summary
if options.showSummary
    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('              PERFORMANCE MODEL SUMMARY             \n');
    fprintf('====================================================\n');
    fprintf('ve (ideal)              : %.6f m/s\n', ve);
    fprintf('c* (ideal)              : %.6f m/s\n', cstar);
    fprintf('Isp (ideal)             : %.6f s\n', Isp);
    fprintf('mdot (ideal)             : %.6f kg/s\n', mdot);
    fprintf('====================================================\n');
    fprintf('\n');

    fprintf('====================================================\n');
    fprintf('                 NOZZLE DESIGN SUMMARY              \n');
    fprintf('====================================================\n');

    fprintf('Type                    : %s\n', type);
    fprintf('CT                      : %.4f [-]\n', CT);
    fprintf('Epsilon                 : %.4f [-]\n', eps);
    fprintf('Lambda                  : %.4f [-]\n', lambdaDiv);

    fprintf('\n');
    fprintf('--------------- Areas ---------------\n');
    fprintf('At                      : %.6e m^2\n', At);
    fprintf('Ae                      : %.6e m^2\n', Ae);
    fprintf('Acc                     : %.6e m^2\n', Acc);

    fprintf('\n');
    fprintf('-------------- Radii ----------------\n');
    fprintf('rt                      : %.6f m\n', rt);
    fprintf('re                      : %.6f m\n', re);
    fprintf('rcc                     : %.6f m\n', rcc);

    fprintf('\n');
    fprintf('------------- Diameters -------------\n');
    fprintf('dt                      : %.6f m\n', 2*rt);
    fprintf('de                      : %.6f m\n', 2*re);
    fprintf('dcc                     : %.6f m\n', 2*rcc);

    fprintf('\n');
    fprintf('------------- Lengths ---------------\n');
    fprintf('Lconv                   : %.6f m\n', lConv);
    fprintf('Ldiv_conical            : %.6f m\n', lDiv);
    fprintf('Ldiv                    : %.6f m\n', lDiv);
    fprintf('Ltotal_nozzle           : %.6f m\n', lTotal);

    fprintf('\n');
    fprintf('-------------- Angles ----------------\n');
    fprintf('alpha                   : %.3f deg\n', alpha_deg);
    fprintf('beta                    : %.3f deg\n', beta_deg);
    fprintf('machCC                     : %.4f [-]\n', machCC);

    fprintf('====================================================\n');
    fprintf('\n');
end
end