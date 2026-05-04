function  [J, massTotal, costTotal] = ccCost(liner, casing, wMass, wCost)
% ccCost - Function to evaluate the weighted cost-mass objective function for the
% combustion chamber sizing optimization.
%
% The combustion chamber is approximated as a cylindrical structure
% consisting of:
%   1) an ablative thermal liner
%   2) a structural metallic casing
%
% The function computes:
%   - liner and casing volumes
%   - corresponding masses
%   - estimated material costs
%   - a normalized weighted objective function used for optimization
%
% INPUT:
%   liner  - struct containing liner properties
%       .thickness     [m]    liner thickness
%       .density       [kg/m^3]
%       .cost          [€/kg]
%
%   casing - struct containing casing properties
%       .thickness     [m]    casing thickness
%       .density       [kg/m^3]
%       .cost          [€/kg]
%
%   wMass  - weighting coefficient for total mass objective [-]
%
%   wCost  - weighting coefficient for total cost objective [-]
%
% OUTPUT:
%   J           - normalized weighted objective function [-]
%
%   massTotal   - total combustion chamber mass [kg]
%
%   costTotal   - total combustion chamber material cost [€]
%
% ASSUMPTIONS:
%   - cylindrical combustion chamber geometry
%   - uniform liner and casing thicknesses
%   - material cost proportional to material mass
%   - only chamber cylindrical section considered
%
% OBJECTIVE FUNCTION:
%
%   J = wMass * (m/m_ref) + wCost * (C/C_ref)
%
% where:
%   m_ref = 100 kg
%   C_ref = 1000 €
%
% EXAMPLE USAGE:
%
%   [J, m, c] = ccCost(liner, casing, 1, 1);

    %% defining constants
    % geometry of grain
    d_grain = 0.89;
    r_grain = d_grain/2;
    l_grain = 1.58;

    % thickness parameters
    t_liner  = liner.thickness;   % [m]
    t_casing = casing.thickness;  % [m]

    % to normalize the cost and mass for optimization
    massRef = 100;     % kg
    costRef = 1000;    % €

    %% getting volumes of each
    v_liner = pi*((r_grain + t_liner)^2 - r_grain^2)*l_grain;
    v_casing = pi*((r_grain + t_liner + t_casing)^2 - (r_grain + t_liner)^2)*l_grain;

    %% masses
    mass_liner = liner.density * v_liner;
    mass_casing = casing.density * v_casing;

    %% defining the cost model
    costTotal = liner.cost * mass_liner + casing.cost * mass_casing;
    massTotal = mass_liner + mass_casing;

    % weigted combination 
    J = wMass * (massTotal / massRef) + wCost * (costTotal / costRef);

end