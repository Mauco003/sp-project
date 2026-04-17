function  [J, massTotal, costTotal] = ccCost(liner, casing, wMass, wCost)
%% -----------------------------------------------------------------------
% Function to peform analysis of cost and mass optimization on the
% combustion chamber sizing
%
% INPUTS:
%
% OUTPUTS:
% 
% UNITS:
%
%
% EXAMPLE USAGE:
%-------------------------------------------------------------------------

%% defining constants
% geometry of grain
d_grain = 0.89;
r_grain = d_grain/2;
l_grain = 1.58;

% thickness parameters
t_liner  = liner.thickness;   % [m]
t_casing = casing.thickness;  % [m]


%% getting volumes of each
v_liner = pi*((r_grain + t_liner)^2 - r_grain^2)*l_grain;
v_casing = pi*((r_grain + t_liner + t_casing)^2 - (r_grain + t_liner)^2)*l_grain;

%% masses
mass_liner = liner.density * v_liner;
mass_casing = casing.density * v_casing;
total_mass = mass_liner + mass_casing;

%% defining the cost model
costTotal = liner.cost * m_liner + casing.cost * m_casing;
massTotal = m_liner + m_casing;

% weigted combination 
J = wMass * massTotal + wCost * costTotal;

end