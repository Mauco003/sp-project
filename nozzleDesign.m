function nozzle = nozzleDesign(At, Ae, rcc, options)
%NOZZLEDESIGN Computes axisymmetric conical nozzle geometry
%
arguments
    At
    Ae
    rcc
    options.alpha = 15 * pi/180
    options.beta  = 30 * pi/180
    options.lambdaDiv = []
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

rt = sqrt(At/pi);
re = sqrt(Ae/pi);

% Chamber area
Acc = pi * rcc^2;

% Conical divergence length
lDiv = (re - rt) / tan(alpha);

% Convergent section length
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
    fprintf('                 NOZZLE DESIGN SUMMARY              \n');
    fprintf('====================================================\n');

    fprintf('Epsilon                 : %.4f [-]\n', Ae/At);
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
    fprintf('Ldiv                    : %.6f m\n', lDiv);
    fprintf('Ltotal_nozzle           : %.6f m\n', lTotal);

    fprintf('\n');
    fprintf('-------------- Angles ----------------\n');
    fprintf('alpha                   : %.3f deg\n', alpha * 180/pi);
    fprintf('beta                    : %.3f deg\n', beta * 180/pi);
    fprintf('====================================================\n\n');
end

end