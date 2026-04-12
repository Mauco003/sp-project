function atm = standardAtmosphere(h)
% STANDARDATMOSPHERE
% International Standard Atmosphere / US Standard Atmosphere style model
% valid up to about 84.852 km
%
% Input:
%   h   - geometric altitude [m]
%
% Output:
%   atm - struct with fields:
%         atm.h     geometric altitude [m]
%         atm.H     geopotential altitude [m]
%         atm.T     temperature [K]
%         atm.p     pressure [Pa]
%         atm.rho   density [kg/m^3]
%         atm.a     speed of sound [m/s]
%
% Supports scalar or vector input.

%% Constants
g0 = 9.80665;          % [m/s^2]
R  = 287.05287;        % [J/(kg*K)]
gamma = 1.4;
Re = 6356766;          % [m] Earth radius for geopotential conversion

%% Layer definitions (geopotential altitude)
Hb = [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852];   % [m]
Lb = [-0.0065, 0, 0.0010, 0.0028, 0, -0.0028, -0.0020];      % [K/m]

% Base temperature and pressure arrays
Tb = zeros(size(Hb));
pb = zeros(size(Hb));

Tb(1) = 288.15;        % [K]
pb(1) = 101325;        % [Pa]

% Precompute base values for each layer
for i = 1:length(Lb)
    if Lb(i) == 0
        Tb(i+1) = Tb(i);
        pb(i+1) = pb(i) * exp(-g0/(R*Tb(i)) * (Hb(i+1)-Hb(i)));
    else
        Tb(i+1) = Tb(i) + Lb(i)*(Hb(i+1)-Hb(i));
        pb(i+1) = pb(i) * (Tb(i+1)/Tb(i))^(-g0/(Lb(i)*R));
    end
end

%% Convert geometric altitude to geopotential altitude
H = Re .* h ./ (Re + h);

%% Allocate outputs
T = zeros(size(H));
p = zeros(size(H));

%% Compute atmosphere point by point
for k = 1:numel(H)

    if H(k) < Hb(1) || H(k) > Hb(end)
        error('Altitude %.2f m is outside model range [0, 84852] m.', h(k));
    end

    % Find layer
    idx = find(H(k) >= Hb, 1, 'last');

    % If exactly at top limit, keep it in previous layer
    if idx == length(Hb)
        idx = length(Hb) - 1;
    end

    H0 = Hb(idx);
    T0 = Tb(idx);
    p0 = pb(idx);
    L  = Lb(idx);

    if L == 0
        T(k) = T0;
        p(k) = p0 * exp(-g0/(R*T0) * (H(k)-H0));
    else
        T(k) = T0 + L*(H(k)-H0);
        p(k) = p0 * (T(k)/T0)^(-g0/(L*R));
    end
end

rho = p ./ (R*T);
a   = sqrt(gamma * R .* T);

%% Output struct
atm.h   = h;
atm.H   = H;
atm.T   = T;
atm.p   = p;
atm.rho = rho;
atm.a   = a;

end