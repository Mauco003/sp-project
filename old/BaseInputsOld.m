function in = BaseInputsOld()
%BASEINPUTS Baseline case for the liquid engine design project
% All units are SI unless otherwise stated.
%
% This file only stores assumptions / baseline data.
% It should NOT perform sizing calculations.
%
% Usage:
%   in = BaseInputs();

%% ------------------------------------------------------------------------
%  CONSTANTS
%  ------------------------------------------------------------------------
in.const.g0  = 9.80665;           % [m/s^2]
in.const.Ru  = 8.314462618e3;     % [J/(kmol*K)]
in.const.atm = 101325;            % [Pa]

%% ------------------------------------------------------------------------
%  REQUIREMENTS / DESIGN POINT
%  ------------------------------------------------------------------------
% may differ
in.req.Thrust      = 100000;                 % [N]
in.req.altitude    = 0;                 % [m] sea level design point
in.req.totalImp    = 25*10e6;       % [Ns]
% for a constant thrust curve
in.req.tburn       = in.req.totalImp/in.req.Thrust;
%in.req.tburn       = 300;                   % [s]

%in.req.maxLossFrac = 0.015;                 % 1.5% max nozzle losses
atm = standardAtmosphere(in.req.altitude);
in.req.pa          = atm.p;
%% ------------------------------------------------------------------------
%  THERMOCHEMICAL DATA
%  ------------------------------------------------------------------------
%  COMPOUNDS

%in.prop.fuel     = 'N2O4';
%in.prop.oxidizer = 'MMH';
%in.prop.name     = 'N2O4_MMH_v2';

%fuel = GetCompoundData(in.prop.fuel);
%ox   = GetCompoundData(in.prop.oxidizer);
%pair = GetPairData(in.prop.name);

% Densities
in.prop.rho_f  = fuel.rho;      % [kg/m^3]
in.prop.rho_ox = ox.rho;        % [kg/m^3]

% PROPELLANT PAIR
% mixture ratio      -> pair.OF_mass
% chamber temperature-> pair.Tc
% molecular weight   -> pair.MW_gas
% gamma              -> pair.k

%in.prop.OF      = pair.OF_mass;   % oxidizer-to-fuel ratio [-]
in.thermo.gamma = pair.k;         % [-]
in.thermo.MW    = pair.MW_gas;    % [kg/kmol]
in.thermo.Tc    = pair.Tc;        % [K]
in.thermo.R     = in.const.Ru / in.thermo.MW;   % [J/(kg*K)]

%% ------------------------------------------------------------------------
%  CHAMBER INPUTS
%  ------------------------------------------------------------------------
% usually this should be given
in.chamber.pc_atm  = 6e6/101325;                       % [atm]
in.chamber.pc      = in.chamber.pc_atm * in.const.atm; % [Pa]

in.injector.fuel.type = 'short_tube_conical';
in.injector.ox.type = 'short_tube_conical';

% desired diameter, define a law into which first guess make
in.injector.fuel.d_mm = 1.5;
in.injector.ox.d_mm = 1.5;

% assume one of the angles
in.injector.alpha_ox_deg = 30;

%% ------------------------------------------------------------------------
%  NOZZLE INPUTS / ASSUMPTIONS
%  ------------------------------------------------------------------------
in.nozzle.type     = 'conical';

% Design exit pressure: sea level optimum in your notes
in.nozzle.pe      = in.req.pa;  % [Pa]
% in.nozzle.eps = 11.2;

%% ------------------------------------------------------------------------
%  PROPELLANT MASS / VOLUME MARGINS
%  ------------------------------------------------------------------------
in.margins.massFrac   = 0.05;   % +5% mass margin
in.margins.volumeFrac = 0.02;   % +2% volume margin

%% ------------------------------------------------------------------------
%  FEED SYSTEM ASSUMPTIONS
%  ------------------------------------------------------------------------
in.feed.type = 'pressure-fed blowdown';

% Assumed propellant bulk velocity in lines
in.feed.u_prop = 10;            % [m/s]

% Assumed pressure losses
in.feed.dp_feed_atm   = 0.5;    % [atm]
in.feed.dp_inj_fracPc = 0.10;   % injector drop = 10% of Pc

%% ------------------------------------------------------------------------
%  TANK ASSUMPTIONS
%  ------------------------------------------------------------------------
in.tanks.geometry    = 'spherical';

% Burst pressure assumption used in the notes
in.tanks.burstFactor = 2.0;     % Pburst = 2 * Ptank
in.tanks.massMargin = 0.05;     % 5%
in.tanks.volumeMargin = 0.02;   % from 2 - 3 %

% Material
in.tanks.material      = 'CarbonFiber';

%% ------------------------------------------------------------------------
% NOZZLE DESIGN CONDITION
% -------------------------------------------------------------------------
% The user may fix either:
%   - exit pressure pe
%   - expansion ratio eps

gamma = in.thermo.gamma;
pc    = in.chamber.pc;

has_pe  = isfield(in.nozzle, 'pe')  && ~isempty(in.nozzle.pe);
has_eps = isfield(in.nozzle, 'eps') && ~isempty(in.nozzle.eps);

if has_pe && ~has_eps
    % pe is fixed -> compute epsilon
    in.nozzle.pe  = in.nozzle.pe;
    in.nozzle.eps = epsilon(gamma, in.nozzle.pe, pc);

elseif ~has_pe && has_eps
    % epsilon is fixed -> compute pe
    in.nozzle.eps = in.nozzle.eps;
    in.nozzle.pe  = peFromEpsilon(in.nozzle.eps, pc, gamma, 'supersonic');

elseif ~has_pe && ~has_eps
    % default: optimum expansion at design point
    in.nozzle.pe  = in.req.pa;
    in.nozzle.eps = epsilon(gamma, in.nozzle.pe, pc);

else
    % both pe and epsilon are given
    epsCheck = epsilon(gamma, in.nozzle.pe, pc);

    if abs(epsCheck - in.nozzle.eps) / in.nozzle.eps > 1e-3
        warning('Given pe and epsilon are not fully consistent with isentropic relation.');
    end
end

end