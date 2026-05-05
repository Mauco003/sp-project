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

%% OPTIONS
wMass = 1;                          % Penalty multiplier on mass
wCost = 1;                          % Penalty multiplier on cost

%% optimization loop
% ranges to search
tL_vals = linspace(1e-5, 50e-2, 1000);   % liner thickness [m]
tC_vals = linspace((pc*grain.dExt0)/(2*hoopStress) * safetyFactor, 0.2, 1000);   % casing thickness [m]

sol.J = 1000000000000;

for i = 1:length(tL_vals)
    for j = 1:length(tC_vals)

        % Set thicknesses
        liner.thickness  = tL_vals(i);
        casing.thickness = tC_vals(j);

        % Run thermal model
        results = combustionChamber(casing, liner, nozzle, propellant, performanceNom, constants, ...
            "makePlot", false, "showSummary", false);

        % Check survival constraint
%        if results.casingSurvives && results.linerSurvives

            % Compute cost + mass
            [J, massTotal, costTotal] = ccCost(liner, casing, wMass, wCost);

            % Store best solution
            if J < sol.J
                sol.J = J;
                sol.tLinerCC = liner.thickness;
                sol.tCasing = casing.thickness;
                sol.mass = massTotal;
                sol.cost = costTotal;
                sol.results = results;
            end
 %       end
    end
end

% liner nozzle thickness
sol.tLinerNozzle = (liner.regressionRate / liner.density * input.burningTime)*2;
end

