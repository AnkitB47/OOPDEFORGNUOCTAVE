function d = computeDistChannel(p, vertices)
    % Define the rectangular channel geometry
    channel_vertices = [0 0; 1 0; 1 -0.5; 2 -0.5; 3 0; 4 0; 4 1; 0 1];

    % Compute distances to the rectangular channel boundary
    d = zeros(size(p, 1), 1);
    for i = 1:size(p, 1)
        d(i) = distToChannel(p(i, :), channel_vertices);
    end
end
