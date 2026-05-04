% Main Water Tank sizing script

% Input parameters
mDot   = 2;                 % [kg/s] mass flow rate of water
pTank  = 10.75e5;           % [Pa] pressure inside the water tank
tc      = 29;                % [s] water supply duration
rho     = 997;               % [kg/m^3] water density
sigma   = 55e6;              % [Pa] allowable material stress

R       = 2077.7;            % [J/(kg*K)] helium gas constant
k       = 1.67;              % [-] helium specific heat ratio

t0He   = 273.15 + 18;       % [K] initial helium temperature
p0He   = 200e5;             % [Pa] initial helium pressure
pFHe   = pTank * 1.5;      % [Pa] final helium pressure

% Water volume calculations

qTot = mDot / rho;         % [m^3/s] volumetric flow rate of water

vW = qTot * tc;            % [m^3] total water volume consumed

vWsf = vW * 1.05;         % [m^3] water volume with 5% safety factor

% Helium thermodynamic calculations

tFHe = t0He * (pFHe / p0He)^((k - 1) / k);          % [K] final helium temperature assuming adiabatic compression/expansion

% Main water tank sizing

% Assume initial helium volume equals 3% of water volume
vTank = vWsf * 1.03;      % [m^3] total tank volume

lTank = 0.5;                % [m] tank length 

rTank = (vTank / (lTank * pi))^(1/2);          % [m]radius of cylinder tank

tTank = pTank * rTank / sigma * 1.5;       % [m] required wall thickness with safety factor 1.5

% Helium tank sizing

vFHe = vWsf * pTank / (p0He - pFHe * (t0He / tFHe));      % [m^3] required final helium volume

m0He = vFHe * p0He / (R * t0He);       % [kg] initial helium mass from ideal gas law

rHe = (3 * vFHe / (4 * pi))^(1/3);        % [m] equivalent spherical helium tank radius

tHe = p0He * rHe / (2 * sigma) * 1.5;     % [m] helium tank wall thickness with safety factor 1.5

%% Results display

fprintf('\n===== MAIN WATER TANK RESULTS =====\n');

fprintf('\nTotal tank volume vTank                = %.6f m^3\n', vTank);

fprintf('Main tank lenght lTank                 = %.6f m\n', lTank);
fprintf('Main tank radius rTank                 = %.6f m\n', rTank);
fprintf('Main tank wall thickness tTank         = %.6f m\n', tTank);

fprintf('\n===== HELIUM TANK RESULTS =====\n');

fprintf('Initial helium mass m0He               = %.6f kg\n', m0He);

fprintf('Helium tank radius rHe                = %.6f m\n', rHe);
fprintf('Helium tank wall thickness tHe        = %.6f m\n', tHe);

fprintf('\n=====================================\n');