function boundaryEdges = calculateBoundaryEdges(convex_hull_points)
    % Calculate the boundary edges from the convex hull points

    % Get the number of points in the convex hull
    numPoints = size(convex_hull_points, 1);

    % Initialize an empty array to store the boundary edges
    boundaryEdges = [];

    % Iterate over each point in the convex hull
    for i = 1:numPoints
        % Get the current point index and the next point index (circular)
        currentIdx = i;
        nextIdx = mod(i, numPoints) + 1;

        % Add the edge formed by connecting the current point to the next point
        boundaryEdges = [boundaryEdges; convex_hull_points(currentIdx, :), convex_hull_points(nextIdx, :)];
    end
end

