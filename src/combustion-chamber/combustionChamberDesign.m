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
% search ranges based on basic thicknesses
results = combustionChamber(casing, liner, nozzle, grain, propellant, performanceNom, constants);

% generating ranges to search based on nominal values            
tL_vals = linspace(results.liner_thickness, 50e-2, 1000);   % liner thickness [m]
tC_vals = linspace((pc*grain.dExt0)/(2*hoopStress) * safetyFactor, 0.2, 1000);   % casing thickness [m]

% artibtrarily large initial value for optimization
sol.J = 10000;

for i = 1:length(tL_vals)
    for j = 1:length(tC_vals)

        % Set thicknesses
        liner.thickness  = tL_vals(i);
        casing.thickness = tC_vals(j);

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
    end
end

% liner nozzle thickness
sol.tLinerNozzle = (liner.regressionRate / liner.density * input.burningTime)*2;

end

