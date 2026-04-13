clearvars
clc


g0 = Constants.g0;

p_c = 70*1e5; % pascal 
T_c = 2342.92; % kelvin
gamma = 1.2386;
a = 1.6177;
n = 0.3816;
rho_ap = 1950;
rho_hptb = 913;
M_M = 22.087;
rho_p = 1/(0.8/rho_ap +  0.2/rho_hptb); % 1859 kg/m3

I_tot = 2.5*1e6;
T = 100*1e3;
t_b = I_tot/T;

p_ext = 101325; % Pa

nozzle = nozzleDesign(T, T_c, p_c, p_ext, gamma, M_M);

m_prop = nozzle.mDot * t_b;

% Pressure [bar]
pressureData = [10.1; 9.7; 10.3; 30.2; 31.0; 29.8; 50.2; 51.0; ...
            50.3; 70.2; 69.0; 69.8; 91.2; 89.1; 89.0];

% Burning rate [mm/s]
burningRateData = [4.0; 3.8; 4.1; 5.6; 6.0; 5.7; 7.0; 7.2; ...
                7.1; 8.4; 8.3; 8.6; 8.8; 9.0; 9.2];

[a, aSigma, n, nSigma, R2] = Uncertainty(pressureData, burningRateData);

r_b = a * (p_c*1e-5)^n;

A_b = nozzle.mDot / (rho_p * r_b*1e-3);

d_i = 203.2; % mm
d_e = 298.45; % mm

L_0 = 1/2 * (3*d_e + d_i); % from Richard Nakka (mm)

% d_p = A_b/(pi * L_0*1e-3)