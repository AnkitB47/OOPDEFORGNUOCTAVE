function generatePolyhedralGeometry(vertices, outputFileName, dimension)
    % Generate polyhedral geometry file based on provided vertices

    % Check if the specified dimension is valid
    if dimension ~= 2 && dimension ~= 3
        error('Invalid dimension. Please specify 2 or 3.');
    end

    % Open the output file for writing
    outputFile = fopen(outputFileName, 'w');
    if outputFile == -1
        error('Error opening the output file for writing.');
    end

    % Write dimension information to the file
    fprintf(outputFile, '# Polyhedral Geometry File\n');
    fprintf(outputFile, '# Dimension (%dd)\n', dimension);
    fprintf(outputFile, 'dimension %d\n\n', dimension);

    % Write vertices to the file
    fprintf(outputFile, '# Vertices (x, y, [z])\n');
    disp(size(vertices));  % Add this line to display the size
    disp(vertices);        % Add this line to display the contents
    for i = 1:size(vertices, 1)
        fprintf(outputFile, 'vertex %f %f', vertices(i, 1), vertices(i, 2));
        if dimension == 3
            fprintf(outputFile, ' %f', vertices(i, 3));
        end
        fprintf(outputFile, '\n');
    end
    fprintf(outputFile, '\n');

    % Calculate edges
    edges = calculateEdges(vertices);

    % Write edges to the file
    fprintf(outputFile, '# Edges (vertex1 vertex2)\n');
    for i = 1:size(edges, 1)
        fprintf(outputFile, 'edge %d %d\n', edges(i, 1), edges(i, 2));
    end
    fprintf(outputFile, '\n');

    % Calculate faces
    faces = calculateFaces(edges, dimension);

    % Write faces to the file
    fprintf(outputFile, '# Faces (vertex1 vertex2 vertex3 ...)\n');
    for i = 1:size(faces, 1)
        fprintf(outputFile, 'face');
        for j = 1:length(faces{i})
            fprintf(outputFile, ' %d', faces{i}(j));
        end
        fprintf(outputFile, '\n');
    end
    fprintf(outputFile, '\n');

    % Close the output file
    fclose(outputFile);

    disp('Polyhedral geometry file generated successfully.');
end

function edges = calculateEdges(vertices)
    % Calculate edges based on vertices

    edges = [];
    n = size(vertices, 1);

    for i = 1:n
        for j = i+1:n
            if isequal(vertices(i, :), vertices(j, :))
                continue;  % Skip self-referencing vertices
            end
            edge = sort([i, j]);
            if ~ismember(edge, edges, 'rows')
                edges = [edges; edge];
            end
        end
    end
end

function faces = calculateFaces(edges, dimension)
    % Calculate faces based on edges

    faces = {};
    while ~isempty(edges)
        currentEdge = edges(1, :);
        connectedEdges = findConnectedEdges(currentEdge, edges);

        if dimension == 2
            face = currentEdge;
        else
            % Concatenate currentEdge and connectedEdges only if dimension is 3
            face = [currentEdge, connectedEdges];
        end

        faces{end+1} = face;  % Use cell array to avoid dimension mismatch
        edges = setdiff(edges, connectedEdges, 'rows');
    end
end

function connectedEdges = findConnectedEdges(edge, edges)
    % Find edges connected to the given edge

    connectedEdges = [];
    for i = 1:size(edges, 1)
        if any(ismember(edge, edges(i, :)))
            connectedEdges = [connectedEdges; edges(i, :)];
        end
    end
end
