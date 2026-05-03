function [P_upstream] = turbopump(q, K_vec, L_vec, D, rho, mu, P_downstream)
% TURBOPUMP compute pressure needed at the upstream of the feeding line
%
% USAGE:
%  [P_upstream] = turbopump(q, K_vec, L_vec, D, rho, mu, P_downstream)
%
% INPUT:
%  q            : Volumetric flow rate (m^3/s)
%  K_vec        : Vector of concentrated loss coefficients (e.g., [K1, K2, K3])
%  L_vec        : Vector of pipe segment lengths (m) (e.g., [L1, L2])
%  D            : Diameter of pipelines (m)
%  rho          : Density of fluid (kg/m^3)
%  mu           : Dynamic viscosity (Pa*s)
%  P_downstream : The pressure at the end of the line (Pa)
% 
% OUTPUT:
%  P_upstram    : Pressure at the start of the line (Pa)

% Calculate Velocity
Area = pi * (D/2)^2;
v = q / Area;

% Reynolds Number
Re = (rho * v * D) / mu;

% Friction Factor (f) - Corrected Logic
if Re <= 2300
    % Laminar Flow
    f = 64 / Re;
else
    % Turbulent Flow (Haaland Equation approximation)
    
    epsilon = 1.5e-6; % assumption 
    f = (-1.8 * log10(((epsilon/D)/3.7)^1.11 + 6.9/Re))^-2;
end

% Distributed Losses (Sum of all pipe segments)
% Total L is the sum of the vector L_vec
L_total = sum(L_vec);
deltaP_distributed = f * (L_total / D) * (rho * v^2 / 2);

% Concentrated Losses (Sum of all K factors)
% Total K is the sum of the vector K_vec
K_total = sum(K_vec);
deltaP_concentrated = K_total * (rho * v^2 / 2);

% Total Pressure Drop and Upstream Pressure
deltaP_tot = deltaP_distributed + deltaP_concentrated;
P_upstream = P_downstream + deltaP_tot;

end