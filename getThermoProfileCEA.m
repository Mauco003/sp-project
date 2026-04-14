function [GasData] = getThermoProfileCEA(wtAP, wtHTPB, p_c_bar, eps_conv, eps_div)
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
%   GasData   - Struct containing arrays of M, T, gamma, cp 
%               aligned from inlet to outlet, plus chamber conditions.
%
% EXAMPLE OF INPUT
%  p_chamber = 70; % bar
%  eps_inlet = 2.0 : -0.01 : 1.01;
%  eps_exit  = 1.01 : 0.01 : 2.0;
%  -------
%  GasInfo = getThermoProfileCEA(80, 20, p_chamber, eps_inlet, eps_exit);
%  -------
%  plot example:
%  plot(GasInfo.Mach, GasInfo.T, 'LineWidth', 2);

    % Ensure eps_conv is descending and eps_div is ascending to avoid CEA crashes
    eps_conv = sort(eps_conv, 'descend');
    eps_div  = sort(eps_div, 'ascend');
    
    N_conv = length(eps_conv);
    N_div  = length(eps_div);
    
    % Preallocation
    T_conv = zeros(1, N_conv); M_conv = zeros(1, N_conv);
    cp_conv = zeros(1, N_conv); gam_conv = zeros(1, N_conv);
    
    T_div = zeros(1, N_div); M_div = zeros(1, N_div);
    cp_div = zeros(1, N_div); gam_div = zeros(1, N_div);

     
    % 1. CONVERGENT SECTION (EQUILIBRIUM)
    
    for i = 1:N_conv
        out = CEA('problem', 'rocket', 'equilibrium', 'frozen', 'nfz', 2, ...
                  'case', 'Stark_Conv', 'p,bar', p_c_bar, 'pi/p', p_c_bar*0.98, ...
                  'o/f', (wtAP/wtHTPB), 'sub', eps_conv(i), ...  
                  'reactants', ...
                      'oxid', 'NH4CLO4(I)', 'wt%', 100, 't,k', 298.15, ...
                      'fuel', 'HTPB', 'C', 7.075, 'H', 10.65, 'O', 0.223, 'N', 0.063, ...
                      'h,J/mol', -58000, 'wt%', 100, 't,k', 298.15, ...
                  'output', 'mks', 'short', 'end');
              
        T_conv(i)   = out.output.eql.temperature(end);
        M_conv(i)   = out.output.eql.mach(end);
        gam_conv(i) = out.output.eql.gamma(end);
        cp_conv(i)  = out.output.eql.cp(end);                  % [kJ/(kg*K)]
    end
    
    % Extract Throat Data from the last run
    T_th   = out.output.eql.temperature(2);
    M_th   = out.output.eql.mach(2);
    gam_th = out.output.eql.gamma(2);
    cp_th  = out.output.eql.cp(2);

    % Extract Chamber Data
    GasData.T_c = out.output.eql.temperature(1);
    GasData.p_c = p_c_bar;

    
    % 2. DIVERGENT SECTION (FROZEN)
    
    for i = 1:N_div
        out = CEA('problem', 'rocket', 'equilibrium', 'frozen', 'nfz', 2, ...
                  'case', 'Stark_Div', 'p,bar', p_c_bar, 'pi/p', p_c_bar*0.98, ...
                  'o/f', (wtAP/wtHTPB), 'sup', eps_div(i), ...  
                  'reactants', ...
                      'oxid', 'NH4CLO4(I)', 'wt%', 100, 't,k', 298.15, ...
                      'fuel', 'HTPB', 'C', 7.075, 'H', 10.65, 'O', 0.223, 'N', 0.063, ...
                      'h,J/mol', -58000, 'wt%', 100, 't,k', 298.15, ...
                  'output', 'mks', 'short', 'end');
              
        T_div(i)   = out.output.froz.temperature(end);
        M_div(i)   = out.output.froz.mach(end);
        gam_div(i) = out.output.froz.gamma(end);
        cp_div(i)  = out.output.froz.cp(end);                  
    end

    
    % 3. DATA ASSEMBLY
    
    % Concatenate arrays to provide a continuous profile from inlet to outlet
    GasData.eps   = [eps_conv, 1.0, eps_div];
    GasData.Mach  = [M_conv, M_th, M_div];
    GasData.T     = [T_conv, T_th, T_div];
    GasData.gamma = [gam_conv, gam_th, gam_div];
    GasData.cp    = [cp_conv, cp_th, cp_div] * 1000;          % Convert kJ to J/(kg*K)
    
end