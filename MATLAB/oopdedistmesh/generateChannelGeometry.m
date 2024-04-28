function [vertices, edges] = generateChannelGeometry(outputFile)
    % Define channel vertices and dimension
    Channel = [
        0, 0;
        1, 0;
        1, -0.5;
        2, -0.5;
        2, 0;
        3, 0;
        3, 1;
        0, 1;
    ];
    dimension = 2;

    % Combine vertices and dimension
    vertices = [Channel, dimension * ones(size(Channel, 1), 1)];

    % Open the output file for writing
    outputFile = fopen(outputFile, 'w');

    if outputFile == -1
        error('Error opening the output file for writing.');
    end

    % Write vertices to the file
    fprintf(outputFile, '# Vertices (x, y, z)\n');
    for i = 1:size(vertices, 1)
        fprintf(outputFile, 'vertex %f %f %f\n', vertices(i, 1), vertices(i, 2), vertices(i, 3));
    end
    fprintf(outputFile, '\n');

    % Define edges
    edges = [
        1, 2;
        2, 3;
        3, 4;
        4, 5;
        5, 6;
        6, 7;
        7, 8;
        8, 1;
    ];

    % Write edges to the file
    fprintf(outputFile, '# Edges (vertex1 vertex2)\n');
    for i = 1:size(edges, 1)
        fprintf(outputFile, 'edge %d %d\n', edges(i, 1), edges(i, 2));
    end
    fprintf(outputFile, '\n');

    % Close the output file
    fclose(outputFile);

    % Plot the geometry
    plotChannelGeometry(vertices, edges);
end

function plotChannelGeometry(vertices, edges)
    figure;
    hold on;

    % Plot edges
    for i = 1:size(edges, 1)
        edge = edges(i, :);
        plot3(vertices(edge, 1), vertices(edge, 2), vertices(edge, 3), 'b-', 'LineWidth', 2);
    end

    % Plot vertices
    scatter3(vertices(:, 1), vertices(:, 2), vertices(:, 3), 'r', 'filled');

    axis equal;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Channel Geometry');

    hold off;
end

