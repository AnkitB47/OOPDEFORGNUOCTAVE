% Generate some random points
points = rand(20,2);

% Find convex hull
[convex_hull_points, stack_indices] = convexHullModified(points);

% Plot original points
scatter(points(:,1), points(:,2), 'filled');
hold on;

% Plot convex hull
plot(convex_hull_points(:,1), convex_hull_points(:,2), 'r-', 'LineWidth', 2);

% Connect outer points directly (last iteration)
plot([convex_hull_points(1,1), convex_hull_points(end,1)], [convex_hull_points(1,2), convex_hull_points(end,2)], 'r-', 'LineWidth', 2);

hold off;
axis equal;

