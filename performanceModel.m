function perf = performanceModel(in)
%PERFORMANCEMODEL Ideal thermodynamic performance model
%
% Current version computes only:
%   - ideal exhaust velocity
%   - ideal characteristic velocity
%   - ideal specific impulse
%
% Inputs required in "in":
%   in.const.g0
%   in.thermo.gamma
%   in.thermo.R
%   in.thermo.Tc
%   in.chamber.pc
%   in.nozzle.pe
%
% Outputs:
%   perf.ve
%   perf.cstar
%   perf.Isp

%% ------------------------------------------------------------------------
% Read inputs
% -------------------------------------------------------------------------
g0    = in.const.g0;
T     = in.req.Thrust;
gamma = in.thermo.gamma;
R     = in.thermo.R;
Tc    = in.thermo.Tc;
pc    = in.chamber.pc;
pe    = in.nozzle.pe;

%% ------------------------------------------------------------------------
% Ideal thermodynamic performance
% -------------------------------------------------------------------------
perf.ve    = exhaustVelocityIdeal(gamma, R, Tc, pe, pc);
perf.cstar = cstarIdeal(R, Tc, gamma);
perf.Isp   = perf.ve / g0;
perf.mdot  = T/(perf.Isp*g0);
%% ------------------------------------------------------------------------
% Store useful extra info
% -------------------------------------------------------------------------
perf.gamma = gamma;
perf.R     = R;
perf.Tc    = Tc;
perf.pc    = pc;
perf.pe    = pe;
 
%% ------------------------------------------------------------------------
% Optional summary
% -------------------------------------------------------------------------
if isfield(in, 'performance') && isfield(in.performance, 'showSummary') ...
        && ~isempty(in.performance.showSummary)
    showSummary = in.performance.showSummary;
else
    showSummary = true;
end

if showSummary
    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('              PERFORMANCE MODEL SUMMARY             \n');
    fprintf('====================================================\n');
    fprintf('ve (ideal)              : %.6f m/s\n', perf.ve);
    fprintf('c* (ideal)              : %.6f m/s\n', perf.cstar);
    fprintf('Isp (ideal)             : %.6f s\n', perf.Isp);
    fprintf('mdot (ideal)             : %.6f kg/s\n', perf.mdot);
    fprintf('====================================================\n');
    fprintf('\n');
end

end