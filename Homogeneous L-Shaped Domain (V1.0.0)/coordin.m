function [boundaryCoords, tangentialUnitVector, normalUnitVector, side] = coordin(longSide_nodes, shortSide_nodes, center, longSide_length, shortSide_length)
% Coordinates calculation of grid points for rectangular conductors

% Inputs:
%   longSide_nodes: Number of nodes on the longer sides
%   shortSide_nodes: Number of nodes on the shorter sides
%   center: Center coordinates of the square domain [x, y]
%   longSide_length: Length of the longer sides
%   shortSide_length: Length of the shorter sides

% Outputs:
%   boundaryCoords: Matrix containing x,y coordinates of nodes (size 2 x total_nodes)
%   normalUnitVector: Vector containing components of normal vectors at points (size 2 x total_nodes)
%   tangentialUnitVector: Vector containing components of tangential vectors at points (size 2 x total_nodes)
%   side: Indicating the side index of each node (ordered from 1 to 4)

%% Section 1: Initialization and Pre-Calculation

% Total number of boundary nodes
total_nodes = 2 * (longSide_nodes + shortSide_nodes);

% Calculate segment lengths for long and short sides
delta_long = longSide_length / longSide_nodes;
delta_short = shortSide_length / shortSide_nodes;

% Check if segment sizes are different
if delta_long ~= delta_short
    disp('Different segment sizes')
end

% Initialize matrix for boundary coordinates and side indices
boundaryCoords = zeros(3, total_nodes);

% Define the four vertices of the rectangle centered at the origin
v1 = [-longSide_length; -shortSide_length] / 2;
v2 = [longSide_length; -shortSide_length] / 2;
v3 = [longSide_length; shortSide_length] / 2;
v4 = [-longSide_length; shortSide_length] / 2;

%% Section 2: Compute Coordinates for Each Side

% Initialize matrices for each side
bottom_side = zeros(3, longSide_nodes);
right_side = zeros(3, shortSide_nodes);
top_side = zeros(3, longSide_nodes);
left_side = zeros(3, shortSide_nodes);

% Bottom side (side 1)
for i = 1:longSide_nodes
    bottom_side(1:2, i) = ((v1 - [delta_long / 2; 0]) + i * [delta_long; 0]);
    bottom_side(3, i) = 1; % side index
end

% Right side (side 2)
for i = 1:shortSide_nodes
    right_side(1:2, i) = ((v2 - [0; delta_short / 2]) + i * [0; delta_short]);
    right_side(3, i) = 2; % side index
end

% Top side (side 3)
for i = 1:longSide_nodes
    top_side(1:2, i) = ((v3 + [delta_long / 2; 0]) - i * [delta_long; 0]);
    top_side(3, i) = 3; % side index
end

% Left side (side 4)
for i = 1:shortSide_nodes
    left_side(1:2, i) = ((v4 + [0; delta_short / 2]) - i * [0; delta_short]);
    left_side(3, i) = 4; % side index
end

% Concatenate all sides into one matrix
local_coords = [bottom_side right_side top_side left_side];

%% Section 3: Calculate Global Coordinates

% Calculate global coordinates of boundary nodes by adding the center offset
boundaryCoords(1:2, :) = local_coords(1:2, :) + kron(center, ones(1, total_nodes));
boundaryCoords(3, :) = local_coords(3, :); % side indices
side = boundaryCoords(3, :); % side vector

%% Section 4: Calculate Unit Vectors

% Initialize unit vectors in global coordinates
tangentialUnitVector = zeros(2, total_nodes);
normalUnitVector = zeros(2, total_nodes);

% Assign unit vectors based on side indices
for i = 1:total_nodes
    if boundaryCoords(3, i) == 1 % Bottom side
        tangentialUnitVector(:, i) = [1; 0]; % Tangential vector along x-axis
        normalUnitVector(:, i) = [0; -1]; % Normal vector pointing down
    elseif boundaryCoords(3, i) == 2 % Right side
        tangentialUnitVector(:, i) = [0; 1]; % Tangential vector along y-axis
        normalUnitVector(:, i) = [1; 0]; % Normal vector pointing right
    elseif boundaryCoords(3, i) == 3 % Top side
        tangentialUnitVector(:, i) = [-1; 0]; % Tangential vector along negative x-axis
        normalUnitVector(:, i) = [0; 1]; % Normal vector pointing up
    elseif boundaryCoords(3, i) == 4 % Left side
        tangentialUnitVector(:, i) = [0; -1]; % Tangential vector along negative y-axis
        normalUnitVector(:, i) = [-1; 0]; % Normal vector pointing left
    end
end

%% Section 5: Extract x, y coordinates
boundaryCoords = boundaryCoords(1:2, :);

%% Section 6: Visualization (Optional)

% % Plot boundary nodes
% figure;
% plot(boundaryCoords(1, :), boundaryCoords(2, :), '*');
% title('Boundary Nodes');
% xlabel('X');
% ylabel('Y');
% % pause;
% 
% % Plot tangential unit vectors
% figure;
% quiver(boundaryCoords(1, :), boundaryCoords(2, :), tangentialUnitVector(1, :), tangentialUnitVector(2, :), 'r');
% xlabel('X');
% ylabel('Y');
% hold on;
% 
% % Plot normal unit vectors
% quiver(boundaryCoords(1, :), boundaryCoords(2, :), normalUnitVector(1, :), normalUnitVector(2, :), 'b');
% legend('Tangential Unit Vectors','Normal Unit Vectors');
% hold off;
% % pause;

end