function plotHollowCylinder(grain)
% plotHollowCylinder - Plots a 3D hollow cylinder with light memory usage.
%
% SYNTAX:
%  plotHollowCylinder(grain)
%
% INPUT:
%  grain - A structure containing the properties of the grain
    
    % Define the number of points for the circle
    numPoints = 36; 
    theta = linspace(0, 2*pi, numPoints);
    
    % Calculate radii
    rExt = grain.dExt0 / 2;
    rInt = grain.dInt0 / 2;
    
    % Matrices for the Outer Cylinder
    xOut = rExt * cos(theta);
    yOut = rExt * sin(theta);
    xOutSurf = [xOut; xOut];
    yOutSurf = [yOut; yOut];
    zOutSurf = [0; grain.L0] * ones(1, numPoints);
    
    % Matrices for the Inner Cylinder
    xInt = rInt * cos(theta);
    yInt = rInt * sin(theta);
    xIntSurf = [xInt; xInt];
    yIntSurf = [yInt; yInt];
    zIntSurf = [0; grain.L0] * ones(1, numPoints);
    
    % Matrices for the top and bottom caps (Annulus)
    rArray = [rInt; rExt];
    xCap = rArray * cos(theta);
    yCap = rArray * sin(theta);
    zCapBottom = zeros(2, numPoints);
    zCapTop = grain.L0 * ones(2, numPoints);
    
    % Plotting
    currentHoldState = ishold;
    hold on;
    
    % Define properties: White faces with colored edges
    surfProps = {'FaceColor', 'w', 'EdgeColor', 'b'};
    
    % Draw the four surfaces
    surf(xOutSurf, yOutSurf, zOutSurf, surfProps{:}); % Outer wall
    surf(xIntSurf, yIntSurf, zIntSurf, surfProps{:}); % Inner wall
    surf(xCap, yCap, zCapBottom, surfProps{:});       % Bottom cap
    surf(xCap, yCap, zCapTop, surfProps{:});          % Top cap
    
    % Formatting
    axis equal;    % Ensures the cylinder doesn't look stretched/oval
    grid on;
    view(3);       % Sets a default 3D isometric view
    
    xlabel('X-Axis');
    ylabel('Y-Axis');
    zlabel('Z-Axis');
    
    % Restore the original hold state
    if ~currentHoldState
        hold off;
    end
end