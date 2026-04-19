clc
clear all

p_chamber = 70; % bar
eps_inlet = 2.0 : -0.1 : 1.01;
eps_exit  = 1.1 : 0.1 : 2.0;

[Chamber, GasInfo] = getThermoProfileCEA_froz(80, 20, p_chamber, eps_inlet, eps_exit);
 
% PROFESSIONAL PLOTTING (T vs eps) - 

% --- Style Settings ---
fontName   = 'Times New Roman';
fontSize   = 12;
techBlue   = [0, 0.4470, 0.7410];       % Deep technical blue
techOrange = [0.8500, 0.3250, 0.0980];  % High-contrast orange
lineColor  = [0.65, 0.65, 0.65];        % Lighter grey for the path
markerSize = 7;                         % Slightly smaller for elegance
lineWidth  = 1.5;                       % Thinner connection lines

% --- Setup Figure ---
f = figure('Color', 'w', 'Position', [100, 100, 800, 500], 'Name', 'Thermal Profile'); 
hold on;

% --- Dynamic Discretization Indices ---
[~, idx_th] = min(GasInfo.eps);
idx_sub = 1 : (idx_th - 1);
idx_sup = (idx_th + 1) : length(GasInfo.eps);

% --- 1. Plot Continuous Path (Grey) ---
plot(GasInfo.eps, GasInfo.T, '-', 'Color', lineColor, ...
     'LineWidth', lineWidth, 'HandleVisibility', 'off');

% --- 2. Plot Subsonic Points (Blue Squares) ---
h_sub = plot(GasInfo.eps(idx_sub), GasInfo.T(idx_sub), 's', ...
             'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', techBlue, ...
             'MarkerSize', markerSize, 'LineWidth', 1, 'DisplayName', 'Subsonic Stations');
         
% --- 3. Plot Supersonic Points (Blue Diamonds) ---
h_sup = plot(GasInfo.eps(idx_sup), GasInfo.T(idx_sup), 'd', ...
             'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', techBlue, ...
             'MarkerSize', markerSize, 'LineWidth', 1, 'DisplayName', 'Supersonic Stations');
         
% --- 4. Plot Throat Point (Orange Circle) ---
h_th = plot(GasInfo.eps(idx_th), GasInfo.T(idx_th), 'o', ...
            'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', techOrange, ...
            'MarkerSize', markerSize + 2, 'LineWidth', 1.2, ...
            'DisplayName', 'Throat (\epsilon = 1)');

% --- Axis Formatting ---
xlabel('Area Ratio, \epsilon [-]', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
ylabel('Static Temperature, T [K]', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');

% Precise Limits with padding
xlim([0.9, max(GasInfo.eps) + 0.1]);

% Dynamic Y-limits with padding to avoid points touching the top/bottom borders
y_min = min(GasInfo.T) - 50;
y_max = max(GasInfo.T) + 50;
ylim([y_min, y_max]);

% Aspect Ratio and Box
box on;     
grid on;

% Refine axes appearance (lighter grid and borders)
ax = gca;
set(ax, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 0.8, ...
         'GridLineStyle', '--', 'GridAlpha', 0.3, 'Layer', 'top');
         
% Legend setup
legend([h_sub, h_th, h_sup], 'Location', 'east', 'Box', 'on', 'FontName', fontName);       
% --- Saving to PDF ---
fileName = 'Nozzle_Thermal_Profile_eps.pdf';
try
    % exportgraphics is the modern standard for cropped, vector output
    exportgraphics(f, fileName, 'ContentType', 'vector');
    fprintf('Figure successfully saved as: %s\n', fileName);
catch
    % Fallback for older MATLAB versions
    print(f, fileName, '-dpdf', '-r300', '-bestfit');
    fprintf('Figure saved using fallback method as: %s\n', fileName);
end