function noz = nozzleDesign(in)
%NOZZLEDESIGN Preliminary nozzle geometry model
%
% This function supports:
%   - conical nozzles
%   - bell nozzles (partially implemented)
%
% Inputs expected in "in":
%   in.req.Thrust
%   in.chamber.pc
%   in.thermo.gamma
%   in.nozzle.CF_opt
%   in.nozzle.epsilon
%   in.nozzle.type
%
% Optional:
%   in.nozzle.alpha_deg
%   in.nozzle.beta_deg
%   in.nozzle.lambdaDiv
%   in.chamber.Mcc
%   in.chamber.r_cc
%   in.nozzle.bellFrac
%   in.nozzle.makePlot
%
% Outputs in "noz":
%   throat/exit areas and radii
%   convergence/divergence lengths
%   lambda
%   chamber area/radius
%   bell preliminary length if requested

%% ------------------------------------------------------------------------
% 0) Read mandatory inputs
% -------------------------------------------------------------------------
Fopt    = in.req.Thrust;
pc      = in.chamber.pc;
gamma   = in.thermo.gamma;
typeStr = lower(in.nozzle.type);

pe = in.nozzle.pe;
epsRatio = in.nozzle.eps;

CFopt = CFideal(gamma, pe, pc);
noz.CT = CFopt;

%% ------------------------------------------------------------------------
% 1) Defaults / optional inputs
% -------------------------------------------------------------------------
if isfield(in.nozzle, 'alpha_deg') && ~isempty(in.nozzle.alpha_deg)

    alpha_deg = in.nozzle.alpha_deg;

elseif isfield(in.nozzle, 'lambdaDiv') && ~isempty(in.nozzle.lambdaDiv)

    alpha_deg = acosd(2*in.nozzle.lambdaDiv - 1);

else

    alpha_deg = 15;   % default assumption

end

if isfield(in.nozzle, 'beta_deg') && ~isempty(in.nozzle.beta_deg)
    beta_deg = in.nozzle.beta_deg;
else
    beta_deg = 30;    % default assumption
end

if isfield(in.chamber, 'Mcc') && ~isempty(in.chamber.Mcc)
    Mcc = in.chamber.Mcc;
else
    Mcc = 0.2;        % default assumption
end

if isfield(in.nozzle, 'bellFrac') && ~isempty(in.nozzle.bellFrac)
    bellFrac = in.nozzle.bellFrac;
else
    bellFrac = 0.60;  % default assumption
end

if isfield(in.nozzle, 'makePlot') && ~isempty(in.nozzle.makePlot)
    makePlot = in.nozzle.makePlot;
else
    makePlot = true;
end
if isfield(in.nozzle, 'showSummary') && ~isempty(in.nozzle.showSummary)
    showSummary = in.nozzle.showSummary;
else
    showSummary = true;
end
%% ------------------------------------------------------------------------
% 2) Throat and exit areas
% -------------------------------------------------------------------------
noz.At = Fopt / (pc * CFopt);
noz.Ae = epsRatio * noz.At;

noz.rt = sqrt(noz.At/pi);
noz.re = sqrt(noz.Ae/pi);

noz.dt = 2*noz.rt;
noz.de = 2*noz.re;

%% ------------------------------------------------------------------------
% 3) Divergence section (conical-equivalent)
% -------------------------------------------------------------------------
noz.alpha_deg = alpha_deg;
noz.alpha_rad = deg2rad(alpha_deg);

if isfield(in.nozzle, 'lambdaDiv') && ~isempty(in.nozzle.lambdaDiv)
    noz.lambdaDiv = in.nozzle.lambdaDiv;
else
    noz.lambdaDiv = (1 + cos(noz.alpha_rad))/2;
end
% Conical divergence length
noz.Ldiv_conical = (noz.re - noz.rt) / tan(noz.alpha_rad);

%% ------------------------------------------------------------------------
% 4) Chamber section from Mcc if r_cc is not given
% -------------------------------------------------------------------------
if isfield(in.chamber, 'r_cc') && ~isempty(in.chamber.r_cc)
    noz.rcc = in.chamber.r_cc;
    noz.Acc = pi * noz.rcc^2;
    noz.Mcc = NaN;   % not used because r_cc was directly imposed
else
    noz.Mcc = Mcc;
    noz.Acc = noz.At * areaRatioIsen(Mcc, gamma);
    noz.rcc = sqrt(noz.Acc/pi);
end

noz.dcc = 2*noz.rcc;

%% ------------------------------------------------------------------------
% 5) Convergent section
% -------------------------------------------------------------------------
noz.beta_deg = beta_deg;
noz.beta_rad = deg2rad(beta_deg);

noz.Lconv = (noz.rcc - noz.rt) / tan(noz.beta_rad);

%% ------------------------------------------------------------------------
% 6) Geometry by nozzle type
% -------------------------------------------------------------------------
switch typeStr

    case 'conical'

        noz.type = 'conical';
        noz.Ldiv = noz.Ldiv_conical;
        noz.Ltotal_nozzle = noz.Lconv + noz.Ldiv;

        % Simple piecewise geometry for plotting
        x_conv = linspace(-noz.Lconv, 0, 100);
        r_conv = linspace(noz.rcc, noz.rt, 100);

        x_div = linspace(0, noz.Ldiv, 150);
        r_div = linspace(noz.rt, noz.re, 150);

        noz.x_profile = [x_conv, x_div];
        noz.r_profile = [r_conv, r_div];

        if makePlot
            figure;
            plot(noz.x_profile,  noz.r_profile, 'LineWidth', 1.5); hold on;
            plot(noz.x_profile, -noz.r_profile, 'LineWidth', 1.5);
            axis equal;
            grid on;
            xlabel('x [m]');
            ylabel('r [m]');
            title('Axisymmetric conical nozzle');
        end

    case 'bell'

        noz.type = 'bell';

        % Bell nozzle length is defined as a fraction of equivalent conical length
        noz.Lbell = bellFrac * noz.Ldiv_conical;
        noz.Ldiv = noz.Lbell;
        noz.Ltotal_nozzle = noz.Lconv + noz.Ldiv;

        % Placeholder fields for future Rao contour implementation
        noz.bellFrac = bellFrac;
        noz.raoImplemented = false;
        noz.x_profile = [];
        noz.r_profile = [];

        % For now: no final contour, only store the target bell length
        % TODO:
        %   - throat blend radius
        %   - Rao initial angle
        %   - Rao exit angle
        %   - parabolic / bell contour generation

        if makePlot
            figure;
            % Plot only convergent section + straight placeholder divergence
            x_conv = linspace(-noz.Lconv, 0, 100);
            r_conv = linspace(noz.rcc, noz.rt, 100);

            x_div = linspace(0, noz.Lbell, 150);
            r_div = linspace(noz.rt, noz.re, 150); % TEMPORARY placeholder only

            x_prof = [x_conv, x_div];
            r_prof = [r_conv, r_div];

            plot(x_prof,  r_prof, 'LineWidth', 1.5); hold on;
            plot(x_prof, -r_prof, 'LineWidth', 1.5);
            axis equal;
            grid on;
            xlabel('x [m]');
            ylabel('r [m]');
            title('Bell nozzle (temporary placeholder contour)');
        end

    otherwise
        error('Unknown nozzle type. Use "conical" or "bell".');
end

%% ------------------------------------------------------------------------
% 7) Summary output
% -------------------------------------------------------------------------
if showSummary

    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('                 NOZZLE DESIGN SUMMARY              \n');
    fprintf('====================================================\n');

    fprintf('Type                    : %s\n', noz.type);
    fprintf('CT                      : %.4f [-]\n', noz.CT);
    fprintf('Epsilon                 : %.4f [-]\n', eps);
    fprintf('Lambda                  : %.4f [-]\n', noz.lambdaDiv);

    fprintf('\n');
    fprintf('--------------- Areas ---------------\n');
    fprintf('At                      : %.6e m^2\n', noz.At);
    fprintf('Ae                      : %.6e m^2\n', noz.Ae);
    fprintf('Acc                     : %.6e m^2\n', noz.Acc);

    fprintf('\n');
    fprintf('-------------- Radii ----------------\n');
    fprintf('rt                      : %.6f m\n', noz.rt);
    fprintf('re                      : %.6f m\n', noz.re);
    fprintf('rcc                     : %.6f m\n', noz.rcc);

    fprintf('\n');
    fprintf('------------- Diameters -------------\n');
    fprintf('dt                      : %.6f m\n', noz.dt);
    fprintf('de                      : %.6f m\n', noz.de);
    fprintf('dcc                     : %.6f m\n', noz.dcc);

    fprintf('\n');
    fprintf('------------- Lengths ---------------\n');
    fprintf('Lconv                   : %.6f m\n', noz.Lconv);
    fprintf('Ldiv_conical            : %.6f m\n', noz.Ldiv_conical);
    fprintf('Ldiv                    : %.6f m\n', noz.Ldiv);
    fprintf('Ltotal_nozzle           : %.6f m\n', noz.Ltotal_nozzle);

    if strcmpi(noz.type, 'bell')
        fprintf('Lbell                   : %.6f m\n', noz.Lbell);
        fprintf('Bell fraction           : %.4f [-]\n', noz.bellFrac);
    end

    fprintf('\n');
    fprintf('-------------- Angles ----------------\n');
    fprintf('alpha                   : %.3f deg\n', noz.alpha_deg);
    fprintf('beta                    : %.3f deg\n', noz.beta_deg);

    if ~isnan(noz.Mcc)
        fprintf('Mcc                     : %.4f [-]\n', noz.Mcc);
    else
        fprintf('Mcc                     : imposed through rcc\n');
    end

    fprintf('====================================================\n');
    fprintf('\n');

end
end