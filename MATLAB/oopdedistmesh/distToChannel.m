function d = distToChannel(p, vertices)
    % Calculate the signed distance from point p to the channel geometry defined by vertices

    num_vertices = size(vertices, 1);
    distances = zeros(num_vertices, 1);

    % Calculate distance to each line segment forming the channel
    for i = 1:num_vertices
        % Get the coordinates of the current and next vertices to form a line segment
        v1 = vertices(i, :);
        v2 = vertices(mod(i, num_vertices) + 1, :);

        % Calculate distance from point p to the line segment defined by v1 and v2
        distances(i) = distToSegment(p, v1, v2);
    end

    % Take the minimum distance as the signed distance to the channel
    d = min(distances);
end
