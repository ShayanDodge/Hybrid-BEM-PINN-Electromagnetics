function [U4, boundarySide4, boundaryNodes] = bemforpinn_func_opt(Neu4, nSide)
% BEM solver for a rectangular subdomain with constant boundary elements.
% Solves the Laplace problem and returns the potential on side 4.
%
% Inputs:
%   Neu4   - Neumann boundary condition on side 4
%            (can be a scalar or a vector)
%   nSide  - Number of boundary elements per side direction
%
% Outputs:
%   U4            - Potential values on side 4
%   boundarySide4 - Coordinates of nodes on side 4
%   boundaryNodes - Coordinates of all boundary nodes

    if nargin == 0
        Neu4 = 0;
        nSide = 15;
    end

    %% Geometry
    center  = [0.75; 0.25];
    lengthX = 0.5;
    lengthY = 0.5;

    nSideX = nSide;
    nSideY = nSide;
    nBoundary = 2 * (nSideX + nSideY);

    elemLength = 2 * (lengthX + lengthY) / nBoundary;

    %% Boundary nodes
    [boundaryNodes, ~, ~, sideLabels] = coordin(nSideX, nSideY, center, lengthX, lengthY);
    boundarySide4 = boundaryNodes(:, sideLabels == 4);

    %% Initialize BEM matrices
    H = zeros(nBoundary, nBoundary);
    G = zeros(nBoundary, nBoundary);

    %% Build H matrix
    parfor i = 1:nBoundary
        for j = 1:nBoundary
            if i == j
                H(i, j) = 0.5;
            else
                H(i, j) = intnormalsquare(sideLabels, i, j, boundaryNodes, center, lengthX, lengthY, []);
            end
        end
    end

    %% Build G matrix
    r1 = elemLength / 2;
    diagCoeff = (r1 / pi) * (log(1 / r1) + 1);

    parfor i = 1:nBoundary
        for j = 1:nBoundary
            if i == j
                G(i, j) = diagCoeff;
            else
                G(i, j) = intfundamsquare(sideLabels, i, j, boundaryNodes, center, lengthX, lengthY, []);
            end
        end
    end

    %% Boundary classification
    indexDir1 = (sideLabels == 1);
    indexDir3 = (sideLabels == 3);
    indexNeu2 = (sideLabels == 2);
    indexNeu4 = (sideLabels == 4);

    ld1 = nnz(indexDir1);
    ld3 = nnz(indexDir3);
    ln2 = nnz(indexNeu2);
    ln4 = nnz(indexNeu4);

    %% Boundary conditions
    U1 = ones(ld1, 1);      % Dirichlet on side 1
    U3 = zeros(ld3, 1);     % Dirichlet on side 3
    Q2 = zeros(ln2, 1);     % Neumann on side 2

    % Neumann condition on side 4
    if isscalar(Neu4)
        Q4 = Neu4 * ones(ln4, 1);
    else
        Q4 = Neu4(:);
        assert(numel(Q4) == ln4, ...
            'Neu4 size mismatch: expected %d values on side 4, but got %d.', ln4, numel(Q4));
    end

    %% Linear system: H*U = G*Q
    A = [H(:, indexNeu2 | indexNeu4), -G(:, indexDir1 | indexDir3)];

    b1 = H(:, indexDir1) * U1;
    b2 = H(:, indexDir3) * U3;
    b3 = G(:, indexNeu2) * Q2;
    b4 = G(:, indexNeu4) * Q4;

    b = -(b1 + b2) + b3 + b4;

    %% Solve system
    Sol = A \ b;

    %% Extract potential on side 4
    U4 = Sol(ln2 + 1 : ln2 + ln4);

end
