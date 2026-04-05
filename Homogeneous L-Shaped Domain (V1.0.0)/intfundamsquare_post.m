function integrale = intfundamsquare_post(sideLabels, sourcePoint, fieldIdx, boundaryNodes, ~, lengthX, lengthY, ~)
% Computes the integral of the fundamental solution over one boundary element
% for post-processing inside a rectangular domain.
%
% Inputs:
%   sideLabels    - Side label for each boundary node
%   sourcePoint   - Observation point inside the domain
%   fieldIdx      - Index of the boundary element
%   boundaryNodes - Coordinates of boundary nodes
%   lengthX       - Rectangle length in x-direction
%   lengthY       - Rectangle length in y-direction
%
% Output:
%   integrale     - Integral of the fundamental solution

    nBoundary = max(size(boundaryNodes));

    %% Integration points on the selected boundary element
    nSub = 20;
    delta = 2 * (lengthX + lengthY) / (nBoundary * nSub);
    points = zeros(2, nSub + 1);

    switch sideLabels(fieldIdx)
        case 1
            startPoint = boundaryNodes(:, fieldIdx) - [nSub/2 * delta; 0];
            for i = 1:nSub + 1
                points(:, i) = startPoint + [(i - 1) * delta; 0];
            end

        case 2
            startPoint = boundaryNodes(:, fieldIdx) - [0; nSub/2 * delta];
            for i = 1:nSub + 1
                points(:, i) = startPoint + [0; (i - 1) * delta];
            end

        case 3
            startPoint = boundaryNodes(:, fieldIdx) + [nSub/2 * delta; 0];
            for i = 1:nSub + 1
                points(:, i) = startPoint - [(i - 1) * delta; 0];
            end

        case 4
            startPoint = boundaryNodes(:, fieldIdx) + [0; nSub/2 * delta];
            for i = 1:nSub + 1
                points(:, i) = startPoint - [0; (i - 1) * delta];
            end
    end

    %% Evaluate integrand
    integrand = zeros(1, nSub + 1);
    sourcePoint = sourcePoint(:);

    for i = 1:nSub + 1
        r2 = sum((sourcePoint - points(:, i)).^2);
        integrand(i) = log(1 / r2) / (4 * pi);
    end

    %% Numerical integration
    integrale = trapz(integrand) * delta;

end
