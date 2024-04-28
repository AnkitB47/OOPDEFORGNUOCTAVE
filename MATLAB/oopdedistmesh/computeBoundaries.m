function boundaries = computeBoundaries(elements)
    % Compute the boundary edges of the mesh given the element connectivity
    % 'elements'.

    % Initialize boundary edge indices
    boundary_edge_indices = [];

    % Iterate over each element
    for i = 1:size(elements, 1)
        % Extract edges of the current element
        edges = [elements(i, [1, 2]); elements(i, [2, 3]); elements(i, [3, 1])];

        % Check each edge
        for j = 1:size(edges, 1)
            edge = edges(j, :);

            % Convert edge to string for comparison
            edge_str = sprintf('%d,%d', edge);

            % Count occurrences of the current edge in the elements matrix
            count = sum(ismember(cellstr(num2str(elements)), edge_str));

            % If the edge occurs only once, it's a boundary edge
            if count == 1
                boundary_edge_indices = [boundary_edge_indices; edge];
            end
        end
    end

    % Remove duplicate boundary edges
    boundary_edge_indices = unique(boundary_edge_indices, 'rows');

    % Store boundary edges in a cell array
    boundaries = {boundary_edge_indices};
end
