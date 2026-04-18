function [V_tank,r_tank] = cooling_tank(Q_tot,t_CJ,SF)
V_tank = Q_tot*t_CJ*SF;
r_tank = (3*V_tank/(4*pi))^(1/3);
end
