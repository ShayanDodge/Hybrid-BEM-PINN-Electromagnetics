function [xGrid, yGrid, potentialField] = bemforpinn_fun_full( ...
    interfaceFluxLeft, interfaceFluxRight, ...
    interfacePotentialBottom, interfacePotentialTop, ...
    nodesPerOuterSide)
%BEM_PINN_SOLVER_FULL Solve the BEM subproblem coupled to a PINN region.
%
%   [xGrid, yGrid, potentialField] = bem_pinn_solver_full( ...
%       interfaceFluxLeft, interfaceFluxRight, ...
%       interfacePotentialBottom, interfacePotentialTop, ...
%       nodesPerOuterSide)
%
% Inputs
% ------
% interfaceFluxLeft
%     Neumann data q = du/dn supplied by the PINN on the left side of
%     the internal BEM-PINN interface.
%
% interfaceFluxRight
%     Neumann data q = du/dn supplied by the PINN on the right side of
%     the internal BEM-PINN interface.
%
% interfacePotentialBottom
%     Dirichlet data u supplied by the PINN on the bottom side of the
%     internal BEM-PINN interface.
%
% interfacePotentialTop
%     Dirichlet data u supplied by the PINN on the top side of the
%     internal BEM-PINN interface.
%
% nodesPerOuterSide
%     Number of constant boundary elements on each outer side.
%     This value must be an even positive integer because each internal
%     interface side uses half as many elements.
%
% Outputs
% -------
% xGrid, yGrid
%     Cartesian grid used to reconstruct the BEM potential field.
%
% potentialField
%     Reconstructed potential in the BEM region. Values inside the PINN
%     subdomain are set to NaN.
%
% Notes
% -----
% The code uses the boundary-integral equation
%
%                   H * u = G * q
%
% where u is the boundary potential and q = du/dn is the outward-normal
% derivative. The geometry and side ordering are inherited from coordin().

%% Validate inputs

validateattributes(nodesPerOuterSide, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, mfilename, 'nodesPerOuterSide');

if mod(nodesPerOuterSide, 2) ~= 0
    error('nodesPerOuterSide must be even.');
end

%% Geometry and boundary discretisation

domainCentre = [0.5, 0.5]';

outerWidth = 1.0;
outerHeight = 1.0;

outerHorizontalElementCount = nodesPerOuterSide;
outerVerticalElementCount = nodesPerOuterSide;

interfaceHorizontalElementCount = outerHorizontalElementCount / 2;
interfaceVerticalElementCount = outerVerticalElementCount / 2;

totalBoundaryElementCount = ...
    2 * (outerHorizontalElementCount + outerVerticalElementCount) + ...
    2 * (interfaceHorizontalElementCount + interfaceVerticalElementCount);

totalBoundaryLength = ...
    2 * (outerWidth + outerHeight) + ...
    2 * (outerWidth / 2 + outerHeight / 2);

boundaryElementLength = ...
    totalBoundaryLength / totalBoundaryElementCount;

[boundaryCoordinates, ~, ~, boundarySideData] = coordin( ...
    outerHorizontalElementCount, ...
    outerVerticalElementCount, ...
    domainCentre, ...
    outerWidth, ...
    outerHeight);

%% Assemble BEM influence matrices

doubleLayerMatrix = zeros( ...
    totalBoundaryElementCount, totalBoundaryElementCount);

singleLayerMatrix = zeros( ...
    totalBoundaryElementCount, totalBoundaryElementCount);

% H matrix: double-layer influence coefficients.
parfor fieldElementIndex = 1:totalBoundaryElementCount
    for sourceElementIndex = 1:totalBoundaryElementCount
        if fieldElementIndex == sourceElementIndex
            doubleLayerMatrix(fieldElementIndex, sourceElementIndex) = 0.5;
        else
            influenceCoefficient = integral_normal_derivative( ...
                boundarySideData, ...
                fieldElementIndex, ...
                sourceElementIndex, ...
                boundaryCoordinates, ...
                domainCentre, ...
                outerWidth, ...
                outerHeight, ...
                []);

            doubleLayerMatrix( ...
                fieldElementIndex, sourceElementIndex) = ...
                influenceCoefficient;
        end
    end
end

% G matrix: single-layer influence coefficients.
halfElementLength = boundaryElementLength / 2;

for fieldElementIndex = 1:totalBoundaryElementCount
    for sourceElementIndex = 1:totalBoundaryElementCount
        if fieldElementIndex == sourceElementIndex
            influenceCoefficient = ...
                (halfElementLength / pi) * ...
                (log(1 / halfElementLength) + 1);
        else
            influenceCoefficient = integral_potential( ...
                boundarySideData, ...
                fieldElementIndex, ...
                sourceElementIndex, ...
                boundaryCoordinates, ...
                domainCentre, ...
                outerWidth, ...
                outerHeight, ...
                []);
        end

        singleLayerMatrix(fieldElementIndex, sourceElementIndex) = ...
            influenceCoefficient;
    end
end

%% Allocate complete boundary fields

boundaryPotential = zeros(totalBoundaryElementCount, 1);
boundaryFlux = zeros(totalBoundaryElementCount, 1);

%% Define boundary-condition masks

% Side ordering inherited from coordin():
%
%     Outer top:       side 3, prescribed potential
%     Outer right:     side 2, prescribed flux
%     Outer bottom:    side 1, prescribed potential
%     Outer left:      side 4, prescribed flux
%
%     Interface bottom: side 5, prescribed potential from PINN
%     Interface right:  side 6, prescribed flux from PINN
%     Interface top:    side 7, prescribed potential from PINN
%     Interface left:   side 8, prescribed flux from PINN

isOuterBottomPotential = false(totalBoundaryElementCount, 1);
isOuterTopPotential = false(totalBoundaryElementCount, 1);
isOuterRightFlux = false(totalBoundaryElementCount, 1);
isOuterLeftFlux = false(totalBoundaryElementCount, 1);

isInterfaceBottomPotential = false(totalBoundaryElementCount, 1);
isInterfaceTopPotential = false(totalBoundaryElementCount, 1);
isInterfaceRightFlux = false(totalBoundaryElementCount, 1);
isInterfaceLeftFlux = false(totalBoundaryElementCount, 1);

outerBoundaryEndIndex = 2 * totalBoundaryElementCount / 3;

isOuterBottomPotential( ...
    1:outerHorizontalElementCount) = true;

isOuterRightFlux( ...
    outerHorizontalElementCount + 1 : ...
    outerHorizontalElementCount + outerVerticalElementCount) = true;

isOuterTopPotential( ...
    outerHorizontalElementCount + outerVerticalElementCount + 1 : ...
    2 * outerHorizontalElementCount + outerVerticalElementCount) = true;

isOuterLeftFlux( ...
    2 * outerHorizontalElementCount + outerVerticalElementCount + 1 : ...
    outerBoundaryEndIndex) = true;

isInterfaceBottomPotential( ...
    outerBoundaryEndIndex + 1 : ...
    outerBoundaryEndIndex + interfaceHorizontalElementCount) = true;

isInterfaceRightFlux( ...
    outerBoundaryEndIndex + interfaceHorizontalElementCount + 1 : ...
    outerBoundaryEndIndex + interfaceHorizontalElementCount + ...
    interfaceVerticalElementCount) = true;

isInterfaceTopPotential( ...
    outerBoundaryEndIndex + interfaceHorizontalElementCount + ...
    interfaceVerticalElementCount + 1 : ...
    outerBoundaryEndIndex + 2 * interfaceHorizontalElementCount + ...
    interfaceVerticalElementCount) = true;

isInterfaceLeftFlux( ...
    outerBoundaryEndIndex + 2 * interfaceHorizontalElementCount + ...
    interfaceVerticalElementCount + 1 : ...
    totalBoundaryElementCount) = true;

%% Prescribed outer-boundary conditions

outerBottomPotential = 1.0;
outerTopPotential = 0.0;
outerRightFlux = 0.0;
outerLeftFlux = 0.0;

%% Convert masks to index vectors

outerBottomPotentialIndices = find(isOuterBottomPotential);
outerTopPotentialIndices = find(isOuterTopPotential);
outerRightFluxIndices = find(isOuterRightFlux);
outerLeftFluxIndices = find(isOuterLeftFlux);

interfaceBottomPotentialIndices = find(isInterfaceBottomPotential);
interfaceTopPotentialIndices = find(isInterfaceTopPotential);
interfaceRightFluxIndices = find(isInterfaceRightFlux);
interfaceLeftFluxIndices = find(isInterfaceLeftFlux);

%% Assign prescribed boundary and interface data

boundaryPotential(outerBottomPotentialIndices) = outerBottomPotential;
boundaryPotential(outerTopPotentialIndices) = outerTopPotential;

boundaryFlux(outerRightFluxIndices) = outerRightFlux;
boundaryFlux(outerLeftFluxIndices) = outerLeftFlux;

boundaryPotential(interfaceBottomPotentialIndices) = ...
    interfacePotentialBottom;

boundaryPotential(interfaceTopPotentialIndices) = ...
    interfacePotentialTop;

boundaryFlux(interfaceRightFluxIndices) = ...
    interfaceFluxRight;

boundaryFlux(interfaceLeftFluxIndices) = ...
    interfaceFluxLeft;

%% Build the linear system for unknown boundary values

knownPotentialIndices = [
    outerBottomPotentialIndices;
    outerTopPotentialIndices;
    interfaceBottomPotentialIndices;
    interfaceTopPotentialIndices
];

knownFluxIndices = [
    outerRightFluxIndices;
    outerLeftFluxIndices;
    interfaceRightFluxIndices;
    interfaceLeftFluxIndices
];

unknownPotentialMatrix = doubleLayerMatrix;
unknownPotentialMatrix(:, knownPotentialIndices) = [];

unknownFluxMatrix = singleLayerMatrix;
unknownFluxMatrix(:, knownFluxIndices) = [];

systemMatrix = [
    unknownPotentialMatrix, ...
    -unknownFluxMatrix
];

knownPotentialContribution = ...
    doubleLayerMatrix(:, knownPotentialIndices) * ...
    boundaryPotential(knownPotentialIndices);

knownFluxContribution = ...
    singleLayerMatrix(:, knownFluxIndices) * ...
    boundaryFlux(knownFluxIndices);

rightHandSide = ...
    -knownPotentialContribution + knownFluxContribution;

%% Solve for unknown boundary potential and flux

boundaryUnknownVector = systemMatrix \ rightHandSide;

%% Reconstruct complete boundary vectors

boundaryPotentialFull = boundaryPotential;
boundaryFluxFull = boundaryFlux;

outerSideElementCount = numel(outerRightFluxIndices);
interfaceSideElementCount = numel(interfaceRightFluxIndices);

% Unknown potential values on Neumann portions.
offset = 0;

boundaryPotentialFull(outerRightFluxIndices) = ...
    boundaryUnknownVector( ...
    offset + 1 : offset + outerSideElementCount);
offset = offset + outerSideElementCount;

boundaryPotentialFull(outerLeftFluxIndices) = ...
    boundaryUnknownVector( ...
    offset + 1 : offset + outerSideElementCount);
offset = offset + outerSideElementCount;

boundaryPotentialFull(interfaceRightFluxIndices) = ...
    boundaryUnknownVector( ...
    offset + 1 : offset + interfaceSideElementCount);
offset = offset + interfaceSideElementCount;

boundaryPotentialFull(interfaceLeftFluxIndices) = ...
    boundaryUnknownVector( ...
    offset + 1 : offset + interfaceSideElementCount);
offset = offset + interfaceSideElementCount;

% Unknown flux values on Dirichlet portions.
boundaryFluxFull(outerBottomPotentialIndices) = ...
    boundaryUnknownVector( ...
    offset + 1 : offset + outerSideElementCount);
offset = offset + outerSideElementCount;

boundaryFluxFull(outerTopPotentialIndices) = ...
    boundaryUnknownVector( ...
    offset + 1 : offset + outerSideElementCount);
offset = offset + outerSideElementCount;

boundaryFluxFull(interfaceBottomPotentialIndices) = ...
    boundaryUnknownVector( ...
    offset + 1 : offset + interfaceSideElementCount);
offset = offset + interfaceSideElementCount;

boundaryFluxFull(interfaceTopPotentialIndices) = ...
    boundaryUnknownVector( ...
    offset + 1 : offset + interfaceSideElementCount);

%% Optional diagnostic plots

showDiagnostics = false;

if showDiagnostics
    figure('Name', 'BEM boundary diagnostics');

    subplot(2, 2, 1);
    plot(boundaryPotentialFull, 'x');
    grid on;
    xlabel('Boundary element index');
    ylabel('u');
    title('Boundary potential');

    subplot(2, 2, 2);
    plot(boundaryFluxFull, 'x');
    grid on;
    xlabel('Boundary element index');
    ylabel('q = \partial u/\partial n');
    title('Boundary flux');

    subplot(2, 2, 3);
    plot(rightHandSide, 'x');
    grid on;
    xlabel('Equation index');
    ylabel('Value');
    title('Linear-system right-hand side');

    subplot(2, 2, 4);
    plot(boundaryUnknownVector, 'x');
    grid on;
    xlabel('Unknown index');
    ylabel('Value');
    title('Solved boundary unknowns');
end

%% Reconstruct the potential field in the BEM region

gridPointCountPerDirection = nodesPerOuterSide;

xCoordinates = linspace( ...
    domainCentre(1) - outerWidth / 2, ...
    domainCentre(1) + outerWidth / 2, ...
    gridPointCountPerDirection);

yCoordinates = linspace( ...
    domainCentre(2) - outerHeight / 2, ...
    domainCentre(2) + outerHeight / 2, ...
    gridPointCountPerDirection);

[xGrid, yGrid] = meshgrid(xCoordinates, yCoordinates);

potentialField = zeros( ...
    gridPointCountPerDirection, gridPointCountPerDirection);

for rowIndex = 2:gridPointCountPerDirection - 1
    for columnIndex = 2:gridPointCountPerDirection - 1
        evaluationPoint = [
            xGrid(rowIndex, columnIndex), ...
            yGrid(rowIndex, columnIndex)
        ];

        singleLayerEvaluation = zeros( ...
            1, totalBoundaryElementCount);

        doubleLayerEvaluation = zeros( ...
            1, totalBoundaryElementCount);

        for boundaryElementIndex = 1:totalBoundaryElementCount
            singleLayerEvaluation(boundaryElementIndex) = ...
                integral_potential_post( ...
                boundarySideData, ...
                evaluationPoint, ...
                boundaryElementIndex, ...
                boundaryCoordinates, ...
                domainCentre, ...
                outerWidth, ...
                outerHeight, ...
                []);

            doubleLayerEvaluation(boundaryElementIndex) = ...
                integral_normal_derivative_post( ...
                boundarySideData, ...
                evaluationPoint, ...
                boundaryElementIndex, ...
                boundaryCoordinates, ...
                domainCentre, ...
                outerWidth, ...
                outerHeight, ...
                []);
        end

        potentialField(rowIndex, columnIndex) = ...
            singleLayerEvaluation * boundaryFluxFull - ...
            doubleLayerEvaluation * boundaryPotentialFull;
    end
end

%% Insert reconstructed outer-boundary values into the field grid

potentialField(:, 1) = ...
    boundaryPotentialFull(flip(outerLeftFluxIndices));

potentialField(:, end) = ...
    boundaryPotentialFull(outerRightFluxIndices);

potentialField(1, 1:end - 1) = outerBottomPotential;

%% Mask the PINN subdomain

interfaceXMinimum = domainCentre(1) - outerWidth / 4;
interfaceXMaximum = domainCentre(1) + outerWidth / 4;
interfaceYMinimum = domainCentre(2) - outerHeight / 4;
interfaceYMaximum = domainCentre(2) + outerHeight / 4;

isInsidePinnDomain = ...
    xGrid > interfaceXMinimum & ...
    xGrid < interfaceXMaximum & ...
    yGrid > interfaceYMinimum & ...
    yGrid < interfaceYMaximum;

potentialField(isInsidePinnDomain) = NaN;

%% Optional potential-field plot

showPotentialField = false;

if showPotentialField
    figure('Name', 'BEM potential field');

    contourf( ...
        xGrid, ...
        yGrid, ...
        potentialField, ...
        200, ...
        'LineColor', ...
        'none');

    axis equal tight;
    colorbar;
    xlabel('x');
    ylabel('y');
    title('BEM potential field');
    colormap turbo;
end

end
