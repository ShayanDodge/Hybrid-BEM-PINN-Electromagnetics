function integrale = intnormalsquare(sideLabels, sourceIdx, fieldIdx, boundaryNodes, ~, lengthX, lengthY, ~)
% Computes the integral of the normal derivative of the fundamental solution
% over one boundary element of a rectangular domain.
%
% Inputs:
%   sideLabels    - Side label for each boundary node
%   sourceIdx     - Index of the source point
%   fieldIdx      - Index of the boundary element
%   boundaryNodes - Coordinates of boundary nodes
%   lengthX       - Rectangle length in x-direction
%   lengthY       - Rectangle length in y-direction
%
% Output:
%   integrale     - Integral of the normal derivative

    nBoundary = max(size(boundaryNodes));
    sourcePoint = boundaryNodes(:, sourceIdx);

    %% Integration points and outward normals
    nSub = 10;
    delta = 2 * (lengthX + lengthY) / (nBoundary * nSub);

    points = zeros(2, nSub + 1);
    normals = zeros(2, nSub + 1);

    switch sideLabels(fieldIdx)
        case 1
            startPoint = boundaryNodes(:, fieldIdx) - [nSub/2 * delta; 0];
            for i = 1:nSub + 1
                points(:, i) = startPoint + [(i - 1) * delta; 0];
                normals(:, i) = [0; -1];
            end

        case 2
            startPoint = boundaryNodes(:, fieldIdx) - [0; nSub/2 * delta];
            for i = 1:nSub + 1
                points(:, i) = startPoint + [0; (i - 1) * delta];
                normals(:, i) = [1; 0];
            end

        case 3
            startPoint = boundaryNodes(:, fieldIdx) + [nSub/2 * delta; 0];
            for i = 1:nSub + 1
                points(:, i) = startPoint - [(i - 1) * delta; 0];
                normals(:, i) = [0; 1];
            end

        case 4
            startPoint = boundaryNodes(:, fieldIdx) + [0; nSub/2 * delta];
            for i = 1:nSub + 1
                points(:, i) = startPoint - [0; (i - 1) * delta];
                normals(:, i) = [-1; 0];
            end
    end

    %% Evaluate integrand
    integrand = zeros(1, nSub + 1);
    sourcePoint = sourcePoint(:);

    for i = 1:nSub + 1
        integrand(i) = normalder(sourcePoint, points(:, i), normals(:, i));
    end

    %% Numerical integration
    integrale = trapz(integrand) * delta;

end