function [exit] = getExitCEA(wtAP, wtHTPB, p_c_bar, eps_exit)
% getExitCEA evaluates the specific flow properties at the nozzle exit 
% using NASA CEA.
% This version assumes FULL FROZEN flow from the combustion chamber (nfz = 1).
%
% INPUTS:
%   wtAP      - Weight percentage of Ammonium Perchlorate [%]
%   wtHTPB    - Weight percentage of HTPB binder [%]
%   p_c_bar   - Chamber pressure [bar]
%   eps_exit  - Area ratio at the nozzle exit (supersonic) [-]
%
% OUTPUTS:
%   exit      - Struct containing flow properties at the nozzle exit 
%               (pressure, mach, cf, cstar).
%
% UNITS:
%   - pressure [Pa]
%   - mach, cf [-]
%   - cstar [m/s]
%  
% EXAMPLE OF INPUT
%  p_chamber = 70; % bar
%  eps_e = 7.67;
%  -------
%  exit_data = getExitCEA(80, 20, p_chamber, eps_e);
%  -------

out = CEA('problem', 'rocket', 'frozen', 'nfz', 1, ...
                  'case', 'Stark_Conv', 'p,bar', p_c_bar, ...
                  'o/f', wtAP/wtHTPB, 'sup', eps_exit, ...  
                  'reactants', ...
                      'oxid', 'NH4CLO4(I)', 'wt%', 100, 't,k', 298.15, ...
                      'fuel', 'HTPB', 'C', 7.075, 'H', 10.65, 'O', 0.223, 'N', 0.063, ...
                      'h,J/mol', -58000, 'wt%', 100, 't,k', 298.15, ...
                  'output', 'mks', 'transport', 'short', 'end');

exit.pressure = out.output.froz.pressure(3)*1e5;
exit.mach = out.output.froz.mach(3);
exit.cf = out.output.froz.cf(3);
exit.cstar = out.output.froz.cstar(3);
end