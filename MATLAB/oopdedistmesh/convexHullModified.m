function [result, stack_indices] = convexHullModified(points)
    n = size(points, 1);
    if n < 3
        error('Convex hull requires at least 3 points');
    end

    % Find the point with the lowest y-coordinate (and lowest x-coordinate if tied)
    [~,idx] = min(points(:,2));
    p = points(idx,:);

    % Sort points based on polar angle in counterclockwise order w.r.t. p
    angles = atan2(points(:,2) - p(2), points(:,1) - p(1));
    [~, order] = sort(angles);
    sorted_points = points(order,:);

    % Sorting via distance
    for i = 1:size(points,1)
      distPoints(i) = norm(points(i,:)-p);
    end

    % Combining the values
    combIndex = [[1:size(points,1)]', angles, distPoints'];

    % Combined Sorting
    sorted_indices = sortrows(combIndex,[2,3]);
    sorted_points = points(sorted_indices(:,1),:);

    % Initialize the stack with the first two sorted points
    stack = sorted_points(1:2,:);

    % Keep track of indices in the stack
    stack_indices = [order(1); order(2)];

    % Graham scan algorithm with modification
    for i = 3:n
        while size(stack, 1) > 1 && ccw(stack(end-1,:), stack(end,:), sorted_points(i,:)) <= 0
            % Check for collinear points and prioritize traversal along them
            while size(stack, 1) > 1 && ccw(stack(end-1,:), stack(end,:), sorted_points(i,:)) == 0
                stack(end,:) = [];
                stack_indices(end) = [];
            end
            stack(end,:) = [];
            stack_indices(end) = [];
        end
        stack = [stack; sorted_points(i,:)];
        stack_indices = [stack_indices; order(i)];
    end

    result = stack;

end

function orientation = ccw(p1, p2, p3)
    orientation = (p2(1) - p1(1)) * (p3(2) - p1(2)) - (p3(1) - p1(1)) * (p2(2) - p1(2));
end

