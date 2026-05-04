function M = machFromAreaRatio(areaRatio, gamma, branch)
% machFromAreaRatio - Computes Mach number from isentropic area ratio A/A*
%
% SYNTAX:
%  M = machFromAreaRatio(areaRatio, gamma, branch)
%
% INPUT:
%  areaRatio - area ratio A/A*  (must be >= 1)
%  gamma     - ratio of specific heats
%  branch    - 'subsonic' or 'supersonic'
%
% OUTPUT:
%  M         - Mach number
%
% Example:
%  Msub = machFromAreaRatio(2.5, 1.4, 'subsonic');
%  Msup = machFromAreaRatio(2.5, 1.4, 'supersonic');

    if areaRatio < 1
        error('areaRatio must be >= 1');
    end

    if nargin < 3
        branch = 'supersonic';
    end

    % Function whose root we want: areaRatioIsen(M,gamma) - areaRatio = 0
    f = @(M) (1 ./ M) .* ...
        ((2/(gamma + 1)) .* (1 + (gamma - 1)/2 .* M.^2)) ...
        .^ ((gamma + 1)/(2*(gamma - 1))) - areaRatio;

    switch lower(branch)
        case 'subsonic'
            % Subsonic solution: 0 < M < 1
            M = fzero(f, [1e-6, 1-1e-6]);

        case 'supersonic'
            % Supersonic solution: M > 1
            M = fzero(f, [1+1e-6, 50]);

        otherwise
            error('branch must be ''subsonic'' or ''supersonic''');
    end

end