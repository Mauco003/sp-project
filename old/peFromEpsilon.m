function pe = peFromEpsilon(epsRatio, pc, gamma, branch)
% PEFROMEPSILON
% Computes exit pressure pe from expansion ratio epsilon, chamber pressure pc,
% and specific heat ratio gamma.
%
% Inputs:
%   epsRatio  - expansion ratio Ae/At [-]
%   pc        - chamber pressure [Pa]
%   gamma     - specific heat ratio [-]
%   branch    - 'supersonic' (default) or 'subsonic'
%
% Output:
%   pe        - exit pressure [Pa]

if nargin < 4 || isempty(branch)
    branch = 'supersonic';
end

if epsRatio <= 1
    error('Expansion ratio epsilon must be > 1.');
end

term1 = ((gamma + 1)/2)^(1/(gamma - 1));

f = @(x) 1 ./ ...
    ( term1 .* x.^(1/gamma) .* ...
    sqrt(((gamma + 1)/(gamma - 1)) .* (1 - x.^((gamma - 1)/gamma))) ) ...
    - epsRatio;

% Search interval in x = pe/pc, with 0 < x < 1
xgrid = [logspace(-8, -1, 2000), linspace(0.1, 0.999999, 3000)];
fgrid = arrayfun(f, xgrid);

% Find sign changes
idx = find(fgrid(1:end-1) .* fgrid(2:end) < 0);

if isempty(idx)
    error('No valid root found for the given epsilon, pc and gamma.');
end

switch lower(branch)
    case 'supersonic'
        k = idx(1);      % first root -> lower pe/pc
    case 'subsonic'
        k = idx(end);    % second root -> higher pe/pc
    otherwise
        error('branch must be "supersonic" or "subsonic".');
end

x1 = xgrid(k);
x2 = xgrid(k+1);

x = fzero(f, [x1, x2]);

pe = x * pc;

end