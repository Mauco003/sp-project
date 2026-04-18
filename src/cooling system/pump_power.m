function [P_elec, P_hyd, P_mech] = pump_power(Q_tot, P_in, P_out, eta_pump, eta_motor,eta_hyd)
    % pump_power: Calculates the electrical power needed for a pump
    %
    % INPUTS:
    %   Q_tot     : Total volumetric flow rate [m^3/s]
    %   P_in      : Inlet pressure [Pa]
    %   P_out     : Required outlet/upstream pressure [Pa]
    %   eta_pump  : Hydraulic/mechanical efficiency of the pump 
    %   eta_motor : Electrical efficiency of the motor 
    %
    % OUTPUTS:
    %   P_elec    : Total electrical power required [W]
    %   P_hyd     : Hydraulic power [W]
    %   P_mech    : Shaft power [W]

    % 1. Calculate the Pressure Differential
    deltaP = P_out - P_in;
    
    if deltaP < 0
        warning('P_out is lower than P_in. No pumping power required.');
        P_hyd = 0; P_mech = 0; P_elec = 0;
        return;
    end

    % 2. Hydraulic Power
    % Equation: P = Q * deltaP
    P_hyd = Q_tot * deltaP/eta_hyd;

    % 3. Mechanical Shaft Power
    % Accounts for internal pump friction and leakage
    P_mech = P_hyd / eta_pump;

    % 4. Total Electrical Power
    % Accounts for electrical losses in the motor/inverter
    P_elec = P_mech / eta_motor;

    % --- Output Report ---
    fprintf('--- PUMP POWER REQUIREMENTS ---\n');
    fprintf('Pressure Increase: %.2f bar\n', deltaP/1e5);
    fprintf('Hydraulic Power:   %.2f kW\n', P_hyd/1000);
    fprintf('Shaft Power:       %.2f kW (at %.0f%% efficiency)\n', P_mech/1000, eta_pump*100);
    fprintf('Electrical Power:  %.2f kW (at %.0f%% motor efficiency)\n\n', P_elec/1000, eta_motor*100);
end