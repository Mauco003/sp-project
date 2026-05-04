function [t, p, rb, grain] = computeBurn(a, n, rhoP, cStar, grain, Athroat)
% computeBurn - Solves the burn ODE for a given set of parameters and returns the burn time, pressure history, burn rate, and updated grain geometry.
% 
% SYNTAX:
%  [t, p, rb, grain] = computeBurn(a, n, rhoP, cStar, grain, Athroat)
%
% INPUT:
%  a        - Vieille's law coefficient [mm/s/bar^n]
%  n        - Vieille's law exponent
%  rhoP     - Propellant density [kg/m^3]
%  cStar    - Characteristic velocity [m/s]
%  grain    - Initial grain geometry struct with fields:
%               .dInt0 - Initial inner diameter [m]
%               .dExt0 - Initial outer diameter [m]
%               .L0 - Initial length [m]
%  Athroat  - Nozzle throat area [m^2]
%
% OUTPUT:
%  t        - Time vector [s]
%  p        - Pressure history in combustion chamber [Pa]
%  rb       - Burn rate history [m/s]
%  grain    - Updated grain geometry struct with fields:
%               .rInt - Inner radius history [m]
%               .L    - Length history [m]


    rExt = grain.dExt0/2;
    rInt0 = grain.dInt0/2;

    x0 = [rInt0; grain.L0];

    SRM.a = a;
    SRM.rhoP = rhoP;
    SRM.cStar = cStar;
    SRM.rExt = rExt;
    SRM.Athroat = Athroat;
    SRM.n = n;

    % Set up ODE options with event function to stop integration when the burn is complete
    options = odeset("Events", @(t, x) eventFunc(t, x, rExt), "RelTol", 1e-9, "AbsTol", 1e-10);

    [t, x] = ode45(@(t, x) burnODE(t, x, SRM), [0, inf], x0, options);

    % Compute burn rate and pressure history from the ODE solution
    grain.rInt = x(:, 1);
    grain.L = x(:, 2);
    Ab = 2*pi*(rExt^2-grain.rInt.^2)+2*pi*grain.L.*grain.rInt;
    p = (a.*rhoP.*cStar.*Ab./Athroat).^(1./(1-n));
    rb = a.*(p.^n);
end

function dx = burnODE(~, x, SRM)
    % burnODE - Defines the system of ODEs for the burn process based on Vieille's law and the geometry of the grain.


    a = SRM.a;
    rhoP = SRM.rhoP;
    cStar = SRM.cStar; % Recovering unit measure balanced to other terms
    rExt = SRM.rExt;
    Athroat = SRM.Athroat;
    n = SRM.n;
    
    rInt = x(1);
    h = x(2);
    Ab = 2*pi*(rExt^2-rInt.^2)+ 2*pi*h.*rInt;
    p = (a.*rhoP.*cStar.*Ab./Athroat).^(1./(1-n));
    rb = a.*(p.^n);

    dx = zeros(2, 1);
    dx(1) = +rb;
    dx(2) = -2*rb;
end

function [value, isterminal, direction] = eventFunc(~, x, rExt)
    % eventFunc - Event function to stop ODE integration when either the inner radius reaches 
    % the outer radius (web burn-through) or the length reaches zero (length burn-through).
    
    value = [rExt - x(1); x(2)]; % Tracks both web burn-through and length burn-through continuously
    isterminal = [1; 1];         % Stop integration if EITHER hits zero
    direction = [0; 0];          % Approach from any direction
end