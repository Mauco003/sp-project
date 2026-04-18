function [ChamberData, GasData] = getThermoProfileCEA_froz(wtAP, wtHTPB, p_c_bar, eps_conv, eps_div)
% getThermoProfileCEA_froz evaluates the thermodynamic and transport 
% properties of the combustion gas along a rocket nozzle using NASA CEA.
% This version assumes FULL FROZEN flow from the combustion chamber (nfz = 1).
%
% INPUTS:
%   wtAP      - Weight percentage of Ammonium Perchlorate [%]
%   wtHTPB    - Weight percentage of HTPB binder [%]
%   p_c_bar   - Chamber pressure [bar]
%   eps_conv  - Array of area ratios for the convergent section (descending)
%   eps_div   - Array of area ratios for the divergent section (ascending)
%
% OUTPUTS:
%   ChamberData - Struct containing thermodynamic properties of the 
%                 combustion chamber (molar_mass, gamma, cp, R_spec, T0, 
%                 conductivity).
%   GasData     - Struct containing arrays of M, T, gamma, cp, 
%                 prandtl, viscosity, conductivity aligned from inlet to outlet.
%
% UNITS:
%   - eps, Mach, gamma, prandtl [-]
%   - T [K]
%   - cp [J/(kg*K)]
%   - viscosity [Pa*s]
%   - conductivity [W/(m*K)]
%   - density of the gas [kg/m^3]
%  
% EXAMPLE OF INPUT
%  p_chamber = 70; % bar
%  eps_inlet = 2.0 : -0.01 : 1.01;
%  eps_exit  = 1.01 : 0.01 : 2.0;
%  -------
%  [ChamberData,GasData] = getThermoProfileCEA_froz(80, 20, p_chamber, eps_inlet, eps_exit);
%  -------
    % Ensure proper sorting to avoid CEA solver issues
    eps_conv = sort(eps_conv, 'descend');
    eps_div  = sort(eps_div, 'ascend');
    
    N_conv = length(eps_conv);
    N_div  = length(eps_div);
    
    % Preallocation
    T_conv = zeros(1, N_conv); M_conv = zeros(1, N_conv);
    cp_conv = zeros(1, N_conv); gam_conv = zeros(1, N_conv);
    pr_conv = zeros(1, N_conv); visc_conv = zeros(1, N_conv);
    cond_conv = zeros(1, N_conv);
    
    T_div = zeros(1, N_div); M_div = zeros(1, N_div);
    cp_div = zeros(1, N_div); gam_div = zeros(1, N_div);
    pr_div = zeros(1, N_div); visc_div = zeros(1, N_div);
    cond_div = zeros(1, N_div);
     
    
    % 1. CONVERGENT SECTION (FULL FROZEN, nfz = 1)
    for i = 1:N_conv
        out = CEA('problem', 'rocket', 'frozen', 'nfz', 1, ...
                  'case', 'Stark_Conv', 'p,bar', p_c_bar, 'pi/p', 69.0846, ...
                  'o/f', (wtAP/wtHTPB), 'sub', eps_conv(i), ...  
                  'reactants', ...
                      'oxid', 'NH4CLO4(I)', 'wt%', 100, 't,k', 298.15, ...
                      'fuel', 'HTPB', 'C', 7.075, 'H', 10.65, 'O', 0.223, 'N', 0.063, ...
                      'h,J/mol', -58000, 'wt%', 100, 't,k', 298.15, ...
                  'output', 'mks', 'transport', 'short', 'end');
              
        T_conv(i)   = out.output.froz.temperature(end);
        M_conv(i)   = out.output.froz.mach(end);
        
        % Extract Gamma and compute Cp via Mayer's Relation
        gam_conv(i) = out.output.froz.gamma(end);
        mw_loc      = out.output.froz.mw(end);
        R_spec      = 8314.46 / mw_loc;
        cp_conv(i)  = (R_spec * gam_conv(i) / (gam_conv(i) - 1)) / 1000; % [kJ/(kg*K)]
        
        % Direct Transport Extraction
        pr_conv(i)   =   out.output.froz.prandtl.froz(end);
        visc_conv(i) =   out.output.froz.viscosity(end);
        cond_conv(i) =   out.output.froz.conduct.froz(end);
    end
    
    % - Extract Chamber Data
    ChamberData.molar_mass = out.output.froz.mw(1);
    ChamberData.R_spec     = 8314.46 / ChamberData.molar_mass;
    ChamberData.gamma      = out.output.froz.gamma(1);
    ChamberData.T          = out.output.froz.temperature(1);
    % Cp via Mayer's relation for Chamber directly in [J/(kg*K)]
    ChamberData.cp         = ChamberData.R_spec * ChamberData.gamma / (ChamberData.gamma - 1); 
    ChamberData.prandtl   =   out.output.froz.prandtl.froz(1);
    ChamberData.viscosity =   out.output.froz.viscosity(1)* 1e-6;
    ChamberData.conductivity = out.output.froz.conduct.froz(1);
    ChamberData.density =  out.output.froz.density(1);

    % - Extract Exit Data
    GasData.exit.cf = out.output.froz.cf(3);
    GasData.exit.cstar = out.output.froz.cstar(3);
    GasData.exit.mach = out.output.froz.mach(3);
    GasData.exit.gamma = out.output.froz.gamma(3);
    GasData.exit.pressure = out.output.froz.pressure(3);
    GasData.exit.temperature =  out.output.froz.temperature(3);

    % - Extract Throat Data
    T_th   = out.output.froz.temperature(2);
    M_th   = out.output.froz.mach(2);
    
    gam_th = out.output.froz.gamma(2);
    mw_th  = out.output.froz.mw(2);
    R_spec_th = 8314.46 / mw_th;
    cp_th  = (R_spec_th * gam_th / (gam_th - 1)) / 1000; % [kJ/(kg*K)]
    
    pr_th   =   out.output.froz.prandtl.froz(2);
    visc_th =   out.output.froz.viscosity(2);
    cond_th =   out.output.froz.conduct.froz(2);
   
    
    % 2. DIVERGENT SECTION (FULL FROZEN, nfz = 1)
    for i = 1:N_div
        out = CEA('problem', 'rocket', 'frozen', 'nfz', 1, ...
                  'case', 'Stark_Div', 'p,bar', p_c_bar, 'pi/p', 69.0846, ...
                  'o/f', (wtAP/wtHTPB), 'sup', eps_div(i), ...  
                  'reactants', ...
                      'oxid', 'NH4CLO4(I)', 'wt%', 100, 't,k', 298.15, ...
                      'fuel', 'HTPB', 'C', 7.075, 'H', 10.65, 'O', 0.223, 'N', 0.063, ...
                      'h,J/mol', -58000, 'wt%', 100, 't,k', 298.15, ...
                  'output', 'mks', 'transport', 'short', 'end');
              
        T_div(i)   = out.output.froz.temperature(end);
        M_div(i)   = out.output.froz.mach(end);
        
        % Extract Gamma and compute Cp via Mayer's Relation
        gam_div(i) = out.output.froz.gamma(end);
        mw_loc     = out.output.froz.mw(end);
        R_spec     = 8314.46 / mw_loc;
        cp_div(i)  = (R_spec * gam_div(i) / (gam_div(i) - 1)) / 1000; % [kJ/(kg*K)]
        
        % Direct Transport Extraction
        pr_div(i)   =   out.output.froz.prandtl.froz(end);
        visc_div(i) =     out.output.froz.viscosity(end);
        cond_div(i) =   out.output.froz.conduct.froz(end);
    end
    
    % 3. DATA ASSEMBLY
    GasData.eps       = [eps_conv, 1.0, eps_div];
    GasData.Mach      = [M_conv, M_th, M_div];
    GasData.T         = [T_conv, T_th, T_div];
    GasData.gamma     = [gam_conv, gam_th, gam_div];
    GasData.cp        = [cp_conv, cp_th, cp_div] * 1000; % Convert kJ to J/(kg*K)
    GasData.prandtl   = [pr_conv, pr_th, pr_div];
    GasData.conductivity = [cond_conv, cond_th, cond_div];
    
   % Cea Matlab normally gives 100 * millipoise unit for viscosity
   % so value is divided for 1e6 to obtain millipoise
       GasData.viscosity = [visc_conv, visc_th, visc_div] * 1e-6; 
    
end