function integral = integral_normal_derivative(side, i, j, ...
    boundaryCoords, ~, longSide_length, shortSide_length, ~)

% integral_normal_derivative calculates the integral using the trapezoidal rule
% to approximate the integral of the normal derivative over a side of a
% polygonal boundary.
%
% Inputs:
%   side: Array indicating the type of side (1, 2, 3, or 4) for boundary conditions.
%   i, j: Indices for boundaryCoords.
%   boundaryCoords: Coordinates of the boundary points.
%   longSide_length, shortSide_length: Lengths of long and short sides.
%
% Outputs:
%   integral: Approximated integral of the normal derivative.

%% Section 1: Initialization

% Calculate the total number of nodes in the boundary
total_nodes = size(boundaryCoords, 2);

% Select the observation point based on the index i
obPoint = boundaryCoords(:, i);

%% Section 2: Determine points for integration and their normals

pr = 4; % Number of points for trapezoidal integration
delta = (2 * (longSide_length + shortSide_length) + ...
    2*(longSide_length/2 + shortSide_length/2)) / total_nodes / pr; % The step size for integration points

% Initialize arrays for integration points and their normals
values = zeros(2, pr + 1); 
normals = zeros(2, pr + 1);

% Determine which side to process based on the index j
id = j; 
switch side(id)
    case 1 % Bottom side
        % Start point for integration
        start = boundaryCoords(:, id) - [pr/2 * delta; 0];
        % Generate integration points and normals
        for k = 1:pr + 1
            values(:, k) = start + [(k-1) * delta; 0];
            normals(:, k) = [0; -1];
        end
    case 2 % Right side
        start = boundaryCoords(:, id) - [0; pr/2 * delta];
        for k = 1:pr + 1
            values(:, k) = start + [0; (k-1) * delta];
            normals(:, k) = [1; 0];
        end
    case 3 % Top side
        start = boundaryCoords(:, id) + [pr/2 * delta; 0];
        for k = 1:pr + 1
            values(:, k) = start - [(k-1) * delta; 0];
            normals(:, k) = [0; 1];
        end
    case 4 % Left side
        start = boundaryCoords(:, id) + [0; pr/2 * delta];
        for k = 1:pr + 1
            values(:, k) = start - [0; (k-1) * delta];
            normals(:, k) = [-1; 0];
        end

    case 5 % Bottom side
        % Start point for integration
        start = boundaryCoords(:, id) - [pr/2 * delta; 0];
        % Generate integration points and normals
        for k = 1:pr + 1
            values(:, k) = start + [(k-1) * delta; 0];
            normals(:, k) = [0; 1];
        end
    case 6 % Right side
        start = boundaryCoords(:, id) - [0; pr/2 * delta];
        for k = 1:pr + 1
            values(:, k) = start + [0; (k-1) * delta];
            normals(:, k) = [-1; 0];
        end
    case 7 % Top side
        start = boundaryCoords(:, id) + [pr/2 * delta; 0];
        for k = 1:pr + 1
            values(:, k) = start - [(k-1) * delta; 0];
            normals(:, k) = [0; -1];
        end
    case 8 % Left side
        start = boundaryCoords(:, id) + [0; pr/2 * delta];
        for k = 1:pr + 1
            values(:, k) = start - [0; (k-1) * delta];
            normals(:, k) = [1; 0];
        end
end

%% Section 3: Optional Visualization

% Uncomment the lines below for visualization purposes
% plot(boundaryCoords(1,:), boundaryCoords(2,:), '*', values(1,:), ...
%     values(2,:), 'o', obPoint(1,1), obPoint(2,1), 'x')
% hold on
% quiver(values(1,:), values(2,:), normals(1,:), normals(2,:))
% hold off

%% Section 4: Calculation of the integral

% Initialize the integrands array
integrands = zeros(1, pr + 1);

% Compute the integrand values using the normal_derivative function
for k = 1:pr + 1
    integrands(k) = normal_derivative(obPoint, values(:, k), normals(:, k));
end

% Approximate the integral using the trapezoidal rule
integral = trapz(integrands) * delta;

end

