% Call the function to generate the rectangular channel geometry
outputFile = 'rectangular_channel_geometry.txt';
generateRectangularChannelGeometry(outputFile);

% Read vertices from the file
vertices = dlmread('rectangular_channel_geometry.txt', ' ', 1, 0);

% Call the function to generate the polyhedral geometry and extract the boundary
[~, ~, convex_hull_points] = generatePolyhedralGeometry(vertices, 'channel_boundary.txt', 2);

% Extract boundary edges
boundaryEdges = calculateBoundaryEdges(convex_hull_points);

% Plot boundary edges
figure;
hold on;
for i = 1:size(boundaryEdges, 1)
    edge = boundaryEdges(i, :);
    plot([convex_hull_points(edge(1), 1), convex_hull_points(edge(2), 1)], ...
         [convex_hull_points(edge(1), 2), convex_hull_points(edge(2), 2)], 'k-', 'LineWidth', 2);
end
axis equal;
xlabel('X');
ylabel('Y');
title('Boundary of Rectangular Channel');
hold off;

disp('Boundary of rectangular channel plotted successfully.');

