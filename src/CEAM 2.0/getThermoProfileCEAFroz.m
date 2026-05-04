function [chamberData, gasData] = getThermoProfileCEAFroz(wtAP, wtHTPB, pcBar, epsConv, epsDiv)
% getThermoProfileCEAFroz - evaluates the thermodynamic and transport 
% properties of the combustion gas along a rocket nozzle using NASA CEA.
% This version assumes FULL FROZEN flow from the combustion chamber (nfz = 1).
%
% SYNTAX:
%  [chamberData, gasData] = getThermoProfileCEAFroz(wtAP, wtHTPB, pcBar, epsConv, epsDiv)
%
% INPUT:
%  wtAP         - Weight percentage of Ammonium Perchlorate [%]
%  wtHTPB       - Weight percentage of HTPB binder [%]
%  pcBar        - Chamber pressure [bar]
%  epsConv      - Array of area ratios for the convergent section (descending)
%  epsDiv       - Array of area ratios for the divergent section (ascending)
%
% OUTPUT:
%  chamberData  - Struct containing thermodynamic properties of the 
%                  combustion chamber (molarMass, gamma, cp, Rspec, T0, 
%                  conductivity).
%  gasData      - Struct containing arrays of M, T, gamma, cp, 
%                  prandtl, viscosity, conductivity aligned from inlet to outlet.
%                  Also exit properties are present (P,T,gamma,cstar, cf, Mach)
%                  EXIT PROPERTIES FOR ADAPTED NOZZLE [!]
% UNIT:
%   - eps, Mach, gamma, prandtl, cf [-]
%   - T [K]
%   - P [Pa]
%   - cp [J/(kg*K)]
%   - viscosity [Pa*s]
%   - conductivity [W/(m*K)]
%   - density of the gas [kg/m^3]
%   - cstar [m/s]
%
%
% EXAMPLE OF INPUT
%  pChamber = 70; % bar
%  epsInlet = 2.0 : -0.01 : 1.01;
%  epsExit  = 1.01 : 0.01 : 2.0;
%  -------
%  [chamberData,gasData] = getThermoProfileCEAFroz(80, 20, pChamber, epsInlet, epsExit);
%  -------
%  NOTE:
% For exit values use gasData.exit.P 

    % Ensure proper sorting to avoid CEA solver issues
    epsConv = sort(epsConv, 'descend');
    epsDiv  = sort(epsDiv, 'ascend');
    
    N_conv = length(epsConv);
    N_div  = length(epsDiv);
    
    % Preallocation
    T_conv = zeros(1, N_conv); M_conv = zeros(1, N_conv);
    cp_conv = zeros(1, N_conv); gam_conv = zeros(1, N_conv);
    pr_conv = zeros(1, N_conv); visc_conv = zeros(1, N_conv);
    cond_conv = zeros(1, N_conv);
    
    T_div = zeros(1, N_div); M_div = zeros(1, N_div);
    cp_div = zeros(1, N_div); gam_div = zeros(1, N_div);
    pr_div = zeros(1, N_div); visc_div = zeros(1, N_div);
    cond_div = zeros(1, N_div);
     
    
    % Convergent section (FULL FROZEN, nfz = 1)
    for i = 1:N_conv
        out = CEA('problem', 'rocket', 'frozen', 'nfz', 1, ...
                  'case', 'Stark_Conv', 'p,bar', pcBar, 'pi/p', 69.0846, ...
                  'o/f', (wtAP/wtHTPB), 'sub', epsConv(i), ...  
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
        Rspec      = 8314.46 / mw_loc;
        cp_conv(i)  = (Rspec * gam_conv(i) / (gam_conv(i) - 1)) / 1000; % [kJ/(kg*K)]
        
        % Direct Transport Extraction
        pr_conv(i)   =   out.output.froz.prandtl.froz(end);
        visc_conv(i) =   out.output.froz.viscosity(end);
        cond_conv(i) =   out.output.froz.conduct.froz(end);
    end
    
    % Extract Chamber Data
    chamberData.molarMass = out.output.froz.mw(1);
    chamberData.Rspec     = 8314.46 / chamberData.molarMass;
    chamberData.gamma      = out.output.froz.gamma(1);
    chamberData.T          = out.output.froz.temperature(1);
    % Cp via Mayer's relation for Chamber directly in [J/(kg*K)]
    chamberData.cp         = chamberData.Rspec * chamberData.gamma / (chamberData.gamma - 1); 
    chamberData.prandtl   =   out.output.froz.prandtl.froz(1);
    chamberData.viscosity =   out.output.froz.viscosity(1)* 1e-6;
    chamberData.conductivity = out.output.froz.conduct.froz(1);
    chamberData.density =  out.output.froz.density(1);

    % Extract Exit Data
    gasData.exit.cf = out.output.froz.cf(3);
    gasData.exit.cstar = out.output.froz.cstar(3);
    gasData.exit.mach = out.output.froz.mach(3);
    gasData.exit.gamma = out.output.froz.gamma(3);
    gasData.exit.pressure = out.output.froz.pressure(3)*1e5;
    gasData.exit.temperature =  out.output.froz.temperature(3);

    % Extract Throat Data
    T_th   = out.output.froz.temperature(2);
    M_th   = out.output.froz.mach(2);
    
    gam_th = out.output.froz.gamma(2);
    mw_th  = out.output.froz.mw(2);
    Rspec_th = 8314.46 / mw_th;
    cp_th  = (Rspec_th * gam_th / (gam_th - 1)) / 1000; % [kJ/(kg*K)]
    
    pr_th   =   out.output.froz.prandtl.froz(2);
    visc_th =   out.output.froz.viscosity(2);
    cond_th =   out.output.froz.conduct.froz(2);
   
    
    % Divergent section (FULL FROZEN, nfz = 1)
    for i = 1:N_div
        out = CEA('problem', 'rocket', 'frozen', 'nfz', 1, ...
                  'case', 'Stark_Div', 'p,bar', pcBar, 'pi/p', 69.0846, ...
                  'o/f', (wtAP/wtHTPB), 'sup', epsDiv(i), ...  
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
        Rspec     = 8314.46 / mw_loc;
        cp_div(i)  = (Rspec * gam_div(i) / (gam_div(i) - 1)) / 1000; % [kJ/(kg*K)]
        
        % Direct Transport Extraction
        pr_div(i)   =   out.output.froz.prandtl.froz(end);
        visc_div(i) =     out.output.froz.viscosity(end);
        cond_div(i) =   out.output.froz.conduct.froz(end);
    end
    
    % Data assembly
    gasData.eps       = [epsConv, 1.0, epsDiv];
    gasData.Mach      = [M_conv, M_th, M_div];
    gasData.T         = [T_conv, T_th, T_div];
    gasData.gamma     = [gam_conv, gam_th, gam_div];
    gasData.cp        = [cp_conv, cp_th, cp_div] * 1000; % Convert kJ to J/(kg*K)
    gasData.prandtl   = [pr_conv, pr_th, pr_div];
    gasData.conductivity = [cond_conv, cond_th, cond_div];
    
   % Cea Matlab normally gives [100 * millipoise] unit for viscosity
   % so value is divided for 1e6 to obtain Poiseuille [Pa * s]
       gasData.viscosity = [visc_conv, visc_th, visc_div] * 1e-6; 
    
end