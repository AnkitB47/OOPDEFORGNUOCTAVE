function [edges, faces] = computeEdgesAndFaces(vertices)
    % Check if the vertices are in 2D or 3D format
    if size(vertices, 2) == 2
        % 2D case
        edges = delaunay(vertices(:, 1), vertices(:, 2));
        faces = [];
    elseif size(vertices, 2) == 3
        % 3D case
        K = convhull(vertices);
        edges = K(:, [1, 2; 2, 3; 3, 1]);
        faces = K;
    else
        % Unsupported dimension
        error('Vertices are not in 2D or 3D format.');
    end

    % Remove duplicate edges
    edges = unique(sort(edges, 2), 'rows');

    % Convert faces to cell array (needed for distmesh2d_polyhedral)
    faces = mat2cell(faces, ones(size(faces, 1), 1), size(faces, 2));
end
