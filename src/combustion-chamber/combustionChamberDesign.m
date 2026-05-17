function [sol] = combustionChamberDesign(input, performanceNom, grain, nozzle, ...
    propellant, casing, liner, constants)

arguments
    input 
    performanceNom 
    grain 
    nozzle 
    propellant 
    casing 
    liner 
    constants Constants = Constants()
end

%% CASING
hoopStress = casing.hoopStress;
safetyFactor = casing.safetyFactor;
pc = input.pcNominal;
casing.thickness = (pc*grain.dExt0)/(2*hoopStress) * safetyFactor;

%% OPTIONS
wMass = 1;                          % Penalty multiplier on mass
wCost = 1;                          % Penalty multiplier on cost

%% optimization loop
% First estimate required liner thickness
initialResults = combustionChamber(casing, liner, nozzle, grain, ...
                                   propellant, performanceNom, constants);

% generating ranges to search based on basic required values
tL_vals = linspace(initialResults.linerRequired, 50e-2, 50);   % [m]
tC_vals = linspace(casing.thickness, 0.2, 50);                 % [m]

% artibtrarily large initial value for optimization
sol.J = 10000;
sol.feasible = false;

for i = 1:length(tL_vals)
    for j = 1:length(tC_vals)

        % Set thicknesses to check
        liner.thickness  = tL_vals(i);
        casing.thickness = tC_vals(j);

        % Recompute thermal performance and survival with new thicknesses
        results = combustionChamber(casing, liner, nozzle, grain, ...
                                    propellant, performanceNom, constants);

        % Only continue if liner and casing survive
        if results.linerSurvives && results.casingSurvives

            % Compute cost + mass based on cost function
            [J, massTotal, costTotal] = ccCost(liner, casing, wMass, wCost);

            % Store best feasible solution
            if J < sol.J
                sol.J = J;
                sol.tLinerCC = liner.thickness;
                sol.tCasing = casing.thickness;
                sol.mass = massTotal;
                sol.cost = costTotal;
                sol.results = results;
                sol.feasible = true;
            end
        end
    end
end


% return if no feasible solution found
if ~sol.feasible
    warning('No feasible combustion chamber design found in the selected search range.');
end

end

