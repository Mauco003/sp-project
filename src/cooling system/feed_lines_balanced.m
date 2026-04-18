function [P_upstream, deltaP_sys, K_orifices, v_tubes, v_final,D,Area_final] = feed_lines_balanced(n_tubes, Q_tot, eps, K_vec, L_vec, rho, mu, P_downstream, nozzle,thick_Af)
    % n_tubes      : Number of tubes
    % Q_tot        : Total flow rate (m^3/s)
    % eps          : Pipe roughness (m)
    % K_vec        : Extra loss coefficients for each tube (bends, valves, etc.)
    % L_vec        : Length of each individual tube (m)
    % rho, mu      : Density and Viscosity of the fluid
    % P_downstream : Pressure at the exit (Pa)
    % Area_final   : Total area where all tubes meet (m^2)

    % Check if we have a length and K value for every tube
    if length(L_vec) ~= n_tubes || length(K_vec) ~= n_tubes
        error('L_vec and K_vec must have the same number of elements as n_tubes.');
    end
    areaEps2 = nozzle.At*2;
    r_Eps2 = sqrt(areaEps2/pi);
    Area_final = ((r_Eps2+thick_Af)^2 - r_Eps2^2 )*pi;
    % 1. Calculate flow and velocity for each tube and the final section
    % We assume each tube gets an equal share of the flow
    Q_tubes = Q_tot / n_tubes;
    Area_tubes = Area_final / n_tubes;
    D = 2 * sqrt(Area_tubes / pi);
    
    v_tubes = Q_tubes / Area_tubes;
    v_final = Q_tot / Area_final; 
    
    % Dynamic pressure inside the tubes
    P_dyn_tubes = 0.5 * rho * v_tubes^2;

    % 2. Calculate Reynolds number
    Re = (rho * v_tubes * D) / mu;

    % 3. Calculate friction factor (f)
    if Re <= 2300
        f = 64 / Re; % Laminar flow
    else
        f = (-1.8 * log10(((eps/D)/3.7)^1.11 + 6.9/Re))^-2; % Turbulent flow
    end

    % 4. Calculate the natural pressure drop for each tube
    deltaP_natural = zeros(1, n_tubes);
    for i = 1:n_tubes
        deltaP_natural(i) = (f * (L_vec(i) / D) + K_vec(i)) * P_dyn_tubes;
    end

    % 5. Find the "hardest" tube (the one with the highest pressure drop)
    % The pump must at least match this one
    [deltaP_max, critical_index] = max(deltaP_natural);
    
    % 6. Balancing: Calculate the extra restriction (Orifice) needed for each tube
    % Short/easy tubes need extra resistance to match the hardest tube
    K_orifices = zeros(1, n_tubes);
    deltaP_missing = zeros(1, n_tubes);
    
    for i = 1:n_tubes
        deltaP_missing(i) = deltaP_max - deltaP_natural(i);
        % Solve for the K value needed to add this missing pressure drop
        K_orifices(i) = deltaP_missing(i) / P_dyn_tubes; 
    end

    % 7. System energy balance
    % The total pressure drop is now equal to the maximum one
    deltaP_sys = deltaP_max;
    
    % Account for the change in velocity at the exit
    P_dyn_change = 0.5 * rho * (v_final^2 - v_tubes^2);
    
    % Calculate the required pump pressure (Upstream)
    P_upstream = (P_downstream + deltaP_sys + P_dyn_change)*1.05; % 

    % --- Simple Output Report ---
    fprintf('--- MANIFOLD BALANCING RESULTS ---\n');
    fprintf('Diametro nominale dei singoli tubi: %.2f mm (%.4f m)\n', D*1000, D);
    fprintf('Critical tube (longest/hardest): Tube %d\n', critical_index);
    fprintf('Balanced system pressure drop: %.2f bar\n', deltaP_sys/1e5);
    fprintf('Required Pump Pressure: %.2f bar\n\n', P_upstream/1e5);
    disp('Extra K (Orifice) to add to each tube:');
    disp(K_orifices);
end