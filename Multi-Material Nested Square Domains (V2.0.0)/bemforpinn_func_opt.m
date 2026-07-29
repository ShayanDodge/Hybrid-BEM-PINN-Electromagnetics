function [ ...
    interfacePotentialRight, ...
    interfacePotentialLeft, ...
    interfacePotentialBottom, ...
    interfacePotentialTop, ...
    interfaceCoordinatesRight, ...
    interfaceCoordinatesLeft, ...
    interfaceCoordinatesBottom, ...
    interfaceCoordinatesTop, ...
    boundaryCoordinates] = ...
    bemforpinn_func_opt( ...
    interfaceFluxRight, ...
    interfaceFluxLeft, ...
    interfaceFluxBottom, ...
    interfaceFluxTop, ...
    nodesPerOuterSide)
%BEM_INTERFACE_POTENTIAL_SOLVER
% Solve the BEM boundary system and recover the potential on the internal
% BEM-PINN interface from prescribed interface fluxes.
%
% Inputs
% ------
% interfaceFluxRight
%     Prescribed normal derivative q = du/dn on the right interface side.
%
% interfaceFluxLeft
%     Prescribed normal derivative q = du/dn on the left interface side.
%
% interfaceFluxBottom
%     Prescribed normal derivative q = du/dn on the bottom interface side.
%
% interfaceFluxTop
%     Prescribed normal derivative q = du/dn on the top interface side.
%
% nodesPerOuterSide
%     Number of constant boundary elements on each outer side. This value
%     must be a positive even integer because each interface side uses half
%     as many elements.
%
% Outputs
% -------
% interfacePotentialRight
%     Computed potential on the right interface side.
%
% interfacePotentialLeft
%     Computed potential on the left interface side.
%
% interfacePotentialBottom
%     Computed potential on the bottom interface side.
%
% interfacePotentialTop
%     Computed potential on the top interface side.
%
% interfaceCoordinatesRight, interfaceCoordinatesLeft,
% interfaceCoordinatesBottom, interfaceCoordinatesTop
%     Coordinates of the boundary-element nodes on each interface side.
%
% boundaryCoordinates
%     Coordinates of all BEM boundary-element nodes.
%
% The boundary integral equation is written as
%
%                       H*u = G*q
%
% where u is the boundary potential and q = du/dn is the outward-normal
% derivative.

%% Validate inputs

validateattributes(nodesPerOuterSide, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'nodesPerOuterSide');

if mod(nodesPerOuterSide, 2) ~= 0
    error('nodesPerOuterSide must be an even integer.');
end

%% Define geometry and boundary discretisation

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

[boundaryCoordinates, ~, ~, boundarySideId] = coordin( ...
    outerHorizontalElementCount, ...
    outerVerticalElementCount, ...
    domainCentre, ...
    outerWidth, ...
    outerHeight);

% Side numbering inherited from coordin():
%   side 5 = bottom interface
%   side 6 = right interface
%   side 7 = top interface
%   side 8 = left interface
interfaceCoordinatesBottom = boundaryCoordinates(:, boundarySideId == 5);
interfaceCoordinatesRight  = boundaryCoordinates(:, boundarySideId == 6);
interfaceCoordinatesTop    = boundaryCoordinates(:, boundarySideId == 7);
interfaceCoordinatesLeft   = boundaryCoordinates(:, boundarySideId == 8);

%% Assemble BEM influence matrices

doubleLayerMatrix = zeros( ...
    totalBoundaryElementCount, totalBoundaryElementCount);

singleLayerMatrix = zeros( ...
    totalBoundaryElementCount, totalBoundaryElementCount);

% Assemble H: double-layer influence coefficients.
parfor fieldElementIndex = 1:totalBoundaryElementCount
    for sourceElementIndex = 1:totalBoundaryElementCount
        if fieldElementIndex == sourceElementIndex
            doubleLayerMatrix( ...
                fieldElementIndex, sourceElementIndex) = 0.5;
        else
            influenceCoefficient = integral_normal_derivative( ...
                boundarySideId, ...
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

% Assemble G: single-layer influence coefficients.
halfElementLength = boundaryElementLength / 2;

parfor fieldElementIndex = 1:totalBoundaryElementCount
    rowValues = zeros(1, totalBoundaryElementCount);

    for sourceElementIndex = 1:totalBoundaryElementCount
        if fieldElementIndex == sourceElementIndex
            influenceCoefficient = ...
                (halfElementLength / pi) * ...
                (log(1 / halfElementLength) + 1);
        else
            influenceCoefficient = integral_potential( ...
                boundarySideId, ...
                fieldElementIndex, ...
                sourceElementIndex, ...
                boundaryCoordinates, ...
                domainCentre, ...
                outerWidth, ...
                outerHeight, ...
                []);
        end

        rowValues(sourceElementIndex) = influenceCoefficient;
    end

    singleLayerMatrix(fieldElementIndex, :) = rowValues;
end

%% Allocate boundary potential and flux vectors

boundaryPotential = zeros(totalBoundaryElementCount, 1);
boundaryFlux = zeros(totalBoundaryElementCount, 1);

%% Define boundary-condition masks

% Outer boundary:
%   side 1 = bottom, prescribed potential
%   side 2 = right,  prescribed flux
%   side 3 = top,    prescribed potential
%   side 4 = left,   prescribed flux
%
% Internal BEM-PINN interface:
%   side 5 = bottom, prescribed flux
%   side 6 = right,  prescribed flux
%   side 7 = top,    prescribed flux
%   side 8 = left,   prescribed flux

isOuterBottomPotential = false(totalBoundaryElementCount, 1);
isOuterTopPotential = false(totalBoundaryElementCount, 1);
isOuterRightFlux = false(totalBoundaryElementCount, 1);
isOuterLeftFlux = false(totalBoundaryElementCount, 1);

isInterfaceBottomFlux = false(totalBoundaryElementCount, 1);
isInterfaceRightFlux = false(totalBoundaryElementCount, 1);
isInterfaceTopFlux = false(totalBoundaryElementCount, 1);
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

isInterfaceBottomFlux( ...
    outerBoundaryEndIndex + 1 : ...
    outerBoundaryEndIndex + interfaceHorizontalElementCount) = true;

isInterfaceRightFlux( ...
    outerBoundaryEndIndex + interfaceHorizontalElementCount + 1 : ...
    outerBoundaryEndIndex + interfaceHorizontalElementCount + ...
    interfaceVerticalElementCount) = true;

isInterfaceTopFlux( ...
    outerBoundaryEndIndex + interfaceHorizontalElementCount + ...
    interfaceVerticalElementCount + 1 : ...
    outerBoundaryEndIndex + 2 * interfaceHorizontalElementCount + ...
    interfaceVerticalElementCount) = true;

isInterfaceLeftFlux( ...
    outerBoundaryEndIndex + 2 * interfaceHorizontalElementCount + ...
    interfaceVerticalElementCount + 1 : ...
    totalBoundaryElementCount) = true;

%% Convert masks to index vectors

outerBottomPotentialIndices = find(isOuterBottomPotential);
outerTopPotentialIndices = find(isOuterTopPotential);
outerRightFluxIndices = find(isOuterRightFlux);
outerLeftFluxIndices = find(isOuterLeftFlux);

interfaceBottomFluxIndices = find(isInterfaceBottomFlux);
interfaceRightFluxIndices = find(isInterfaceRightFlux);
interfaceTopFluxIndices = find(isInterfaceTopFlux);
interfaceLeftFluxIndices = find(isInterfaceLeftFlux);

%% Assign prescribed boundary data

outerBottomPotential = 1.0;
outerTopPotential = 0.0;
outerRightFlux = 0.0;
outerLeftFlux = 0.0;

boundaryPotential(outerBottomPotentialIndices) = outerBottomPotential;
boundaryPotential(outerTopPotentialIndices) = outerTopPotential;

boundaryFlux(outerRightFluxIndices) = outerRightFlux;
boundaryFlux(outerLeftFluxIndices) = outerLeftFlux;

boundaryFlux(interfaceBottomFluxIndices) = interfaceFluxBottom;
boundaryFlux(interfaceRightFluxIndices) = interfaceFluxRight;
boundaryFlux(interfaceTopFluxIndices) = interfaceFluxTop;
boundaryFlux(interfaceLeftFluxIndices) = interfaceFluxLeft;

%% Build the linear system for unknown boundary values

knownPotentialIndices = [
    outerBottomPotentialIndices;
    outerTopPotentialIndices
];

knownFluxIndices = [
    outerRightFluxIndices;
    outerLeftFluxIndices;
    interfaceRightFluxIndices;
    interfaceLeftFluxIndices;
    interfaceBottomFluxIndices;
    interfaceTopFluxIndices
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

%% Solve for unknown boundary values

boundaryUnknownVector = systemMatrix \ rightHandSide;

%% Extract interface potentials from the solved vector

% Because only sides 1 and 3 have prescribed potential, the first block of
% boundaryUnknownVector contains the unknown potential on all Neumann sides
% in the same order in which those sides occur after deleting columns 1 and
% 3 from H:
%
%   outer right, outer left, interface bottom, interface right,
%   interface top, interface left.
%
% The original implementation extracted the interface values by fixed
% offsets. Here the offsets are computed explicitly from the side sizes.

outerRightCount = numel(outerRightFluxIndices);
outerLeftCount = numel(outerLeftFluxIndices);
interfaceBottomCount = numel(interfaceBottomFluxIndices);
interfaceRightCount = numel(interfaceRightFluxIndices);
interfaceTopCount = numel(interfaceTopFluxIndices);
interfaceLeftCount = numel(interfaceLeftFluxIndices);

offset = outerRightCount + outerLeftCount;

interfacePotentialBottom = boundaryUnknownVector( ...
    offset + 1 : offset + interfaceBottomCount);
offset = offset + interfaceBottomCount;

interfacePotentialRight = boundaryUnknownVector( ...
    offset + 1 : offset + interfaceRightCount);
offset = offset + interfaceRightCount;

interfacePotentialTop = boundaryUnknownVector( ...
    offset + 1 : offset + interfaceTopCount);
offset = offset + interfaceTopCount;

interfacePotentialLeft = boundaryUnknownVector( ...
    offset + 1 : offset + interfaceLeftCount);

end