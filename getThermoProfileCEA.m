function [ChamberData,GasData] = getThermoProfileCEA(wtAP, wtHTPB, p_c_bar, eps_conv, eps_div)
% getThermoProfileCEA evaluates the thermodynamic properties of the 
% combustion gas along a rocket nozzle using NASA CEA.
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
%                 combustion chamber (molar_mass, gamma, cp, R_spec).
%   GasData     - Struct containing arrays of M, T, gamma, cp, 
%                 prandtl, and viscosity aligned from inlet to outlet.
%
%  - Units
%    - eps, Mach, gamma, prandtl [-]
%    - T [K]
%    - cp [J/(kg*K)]
%    - viscosity [Pa*s]
%     
%
% EXAMPLE OF INPUT
%  p_chamber = 70; % bar
%  eps_inlet = 2.0 : -0.01 : 1.01;
%  eps_exit  = 1.01 : 0.01 : 2.0;
%  -------
%  [ChamberData,GasData] = getThermoProfileCEA(80, 20, p_chamber, eps_inlet, eps_exit);
%  -------


    % Ensure eps_conv is descending and eps_div is ascending to avoid CEA crashes
    eps_conv = sort(eps_conv, 'descend');
    eps_div  = sort(eps_div, 'ascend');
    
    N_conv = length(eps_conv);
    N_div  = length(eps_div);
    
    % Preallocation
    T_conv = zeros(1, N_conv); M_conv = zeros(1, N_conv);
    cp_conv = zeros(1, N_conv); gam_conv = zeros(1, N_conv);
    pr_conv = zeros(1, N_conv); visc_conv = zeros(1, N_conv);
    
    T_div = zeros(1, N_div); M_div = zeros(1, N_div);
    cp_div = zeros(1, N_div); gam_div = zeros(1, N_div);
    pr_div = zeros(1, N_div); visc_div = zeros(1, N_div);
     
    % 1. CONVERGENT SECTION (EQUILIBRIUM + ANALYTICAL FROZEN EXTRACTION)
    for i = 1:N_conv
        out = CEA('problem', 'rocket', 'equilibrium', 'frozen', 'nfz', 2, ...
                  'case', 'Stark_Conv', 'p,bar', p_c_bar, 'pi/p', 69.0846, ...
                  'o/f', (wtAP/wtHTPB), 'sub', eps_conv(i), ...  
                  'reactants', ...
                      'oxid', 'NH4CLO4(I)', 'wt%', 100, 't,k', 298.15, ...
                      'fuel', 'HTPB', 'C', 7.075, 'H', 10.65, 'O', 0.223, 'N', 0.063, ...
                      'h,J/mol', -58000, 'wt%', 100, 't,k', 298.15, ...
                  'output', 'mks', 'transport', 'end');
              
        T_conv(i)   = out.output.eql.temperature(end);
        M_conv(i)   = out.output.eql.mach(end);
        
        
        % Analytical extraction of Thermo Gamma and Frozen Cp via Mayer's relation
        gam_s_loc   = out.output.eql.gamma(end);
        dlvpt_loc   = out.output.eql.dlvpt(end);
        mw_loc      = out.output.eql.mw(end);
        R_spec      = 8314.46 / mw_loc;
        
        gam_conv(i) = - gam_s_loc * dlvpt_loc;
        cp_conv(i)  = (R_spec * gam_conv(i) / (gam_conv(i) - 1)) / 1000; % [kJ/(kg*K)]
        
        pr_conv(i)   = out.output.eql.prandtl.froz(end);
        visc_conv(i) = out.output.eql.viscosity(end);
    end
    % Extract Chamber Data 
    R_spec_ch = 8314.46 / out.output.eql.mw(1);
    ChamberData.molar_mass = out.output.eql.mw(1);
    ChamberData.gamma = - out.output.eql.gamma(1) * out.output.eql.dlvpt(1);
    ChamberData.cp  = (R_spec_ch * ChamberData.gamma / (ChamberData.gamma - 1)) / 1000;
    ChamberData.R_spec = R_spec;


    % Extract Throat Data (Analytical Extraction)
    T_th      = out.output.eql.temperature(2);
    M_th      = out.output.eql.mach(2);
    
    R_spec_th = 8314.46 / out.output.eql.mw(2);
    gam_th    = - out.output.eql.gamma(2) * out.output.eql.dlvpt(2);
    cp_th     = (R_spec_th * gam_th / (gam_th - 1)) / 1000;
    
    pr_th     = out.output.eql.prandtl.froz(2);
    visc_th   = out.output.eql.viscosity(2);
 
   
    % 2. DIVERGENT SECTION (FROZEN)
    for i = 1:N_div
        out = CEA('problem', 'rocket', 'equilibrium', 'frozen', 'nfz', 2, ...
                  'case', 'Stark_Div', 'p,bar', p_c_bar, 'pi/p', 69.0846, ...
                  'o/f', (wtAP/wtHTPB), 'sup', eps_div(i), ...  
                  'reactants', ...
                      'oxid', 'NH4CLO4(I)', 'wt%', 100, 't,k', 298.15, ...
                      'fuel', 'HTPB', 'C', 7.075, 'H', 10.65, 'O', 0.223, 'N', 0.063, ...
                      'h,J/mol', -58000, 'wt%', 100, 't,k', 298.15, ...
                  'output', 'mks', 'transport', 'end');
              
        T_div(i)   = out.output.froz.temperature(end);
        M_div(i)   = out.output.froz.mach(end);
        gam_div(i) = out.output.froz.gamma(end);
        cp_div(i)  = out.output.froz.cp(end);    
        
        pr_div(i)   = out.output.eql.prandtl.froz(end);
        visc_div(i) = out.output.eql.viscosity(end);                              
    end
    
    % 3. DATA ASSEMBLY
    % Concatenate arrays to provide a continuous profile from inlet to outlet
    GasData.eps       = [eps_conv, 1.0, eps_div];
    GasData.Mach      = [M_conv, M_th, M_div];
    GasData.T         = [T_conv, T_th, T_div];
    GasData.gamma     = [gam_conv, gam_th, gam_div];
    GasData.cp        = [cp_conv, cp_th, cp_div] * 1000; % Convert kJ to J/(kg*K)
    GasData.prandtl   = [pr_conv, pr_th, pr_div];
    GasData.viscosity = [visc_conv, visc_th, visc_div] * 1e-6;
    % NOTE
    % Cea Matlab normally gives 100 * millipoise unit for viscosity
    % so value is divided for 1e6 to obtain millipoise
end