function [inflowIndices, outflowIndices, boundaryIndices] = convexHullMesh(points)
    % Initialize arrays for inflow, outflow, and boundary indices
    inflowIndices = [];
    outflowIndices = [];
    boundaryIndices = [];

    % Find the convex hull of the points
    K = convhull(points);

    % Iterate through each edge of the convex hull
    for i = 1:size(K, 1)
        % Get the indices of the edge vertices
        idx1 = K(i, 1);
        idx2 = K(i, 2);

        % Check if the edge crosses the cavities
        if ~(isInCavity(points(idx1, :)) && isInCavity(points(idx2, :)))
            % If neither point is in the cavity, it's part of the boundary
            boundaryIndices = [boundaryIndices; idx1; idx2];
        else
            % If one point is in the cavity and the other is not, it's an inflow/outflow point
            if isInCavity(points(idx1, :))
                inflowIndices = [inflowIndices; idx1];
                outflowIndices = [outflowIndices; idx2];
            else
                inflowIndices = [inflowIndices; idx2];
                outflowIndices = [outflowIndices; idx1];
            end
        end
    end
end

function inside = isInCavity(point)
    % Check if the point is inside either of the rectangular cavities
    inside = (point(1) > -0.55 && point(1) < -0.45 && point(2) > -1 && point(2) < 0) || ...
             (point(1) > 0.45 && point(1) < 0.55 && point(2) > -1 && point(2) < 0);
end

