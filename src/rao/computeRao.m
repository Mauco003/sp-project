function [outputArg1,outputArg2] = computeRao(At, epsilon, k, theta)
% computeRao - This function computes the parameters for a RAO nozzle
%
% Input Arguments
%   At    [m^2] - Required throat area
%   epsilon [-] - Expansion ratio
%   k       [-] - Efficiency fraction compared to a 15deg conical nozzle
%   theta [rad] - Vector containing [convergent exit] angles

Ae = At*epsilon;
rt = sqrt(At/pi);
re = sqrt(Ae/pi);


l = k*(re-rt)/tan(15*pi/180);

% Parabolic section



end

