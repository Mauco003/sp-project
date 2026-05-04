function [exit] = getExitCEA(wtAP, wtHTPB, pcBar, epsExit)
% getExitCEA  evaluates the specific flow properties at the nozzle exit 
% using NASA CEA.
% This version assumes FULL FROZEN flow from the combustion chamber (nfz = 1).
%
% SYNTAX:
%  [exit] = getExitCEA(wtAP, wtHTPB, pcBar, epsExit)
%
% INPUT:
%  wtAP      - Weight percentage of Ammonium Perchlorate [%]
%  wtHTPB    - Weight percentage of HTPB binder [%]
%  pcBar   - Chamber pressure [bar]
%  epsExit  - Area ratio at the nozzle exit (supersonic) [-]
%
% OUTPUT:
%  exit      - Struct containing flow properties at the nozzle exit 
%               (pressure, mach, cf, cstar).
%
% UNIT:
%  - pressure [Pa]
%  - mach, cf [-]
%  - cstar [m/s]
%  
% EXAMPLE OF INPUT
%  pChamber = 70; % bar
%  epsE = 7.67;
%  -------
%  exit_data = getExitCEA(80, 20, pChamber, epsE);
%  -------

    % CEA call
    out = CEA('problem', 'rocket', 'frozen', 'nfz', 1, ...
                    'case', 'Stark_Conv', 'p,bar', pcBar, ...
                    'o/f', wtAP/wtHTPB, 'sup', epsExit, ...  
                    'reactants', ...
                        'oxid', 'NH4CLO4(I)', 'wt%', 100, 't,k', 298.15, ...
                        'fuel', 'HTPB', 'C', 7.075, 'H', 10.65, 'O', 0.223, 'N', 0.063, ...
                        'h,J/mol', -58000, 'wt%', 100, 't,k', 298.15, ...
                    'output', 'mks', 'transport', 'short', 'end');

    % Extract exit properties
    exit.pressure = out.output.froz.pressure(3)*1e5;
    exit.mach = out.output.froz.mach(3);
    exit.cf = out.output.froz.cf(3);
    exit.cstar = out.output.froz.cstar(3);
end