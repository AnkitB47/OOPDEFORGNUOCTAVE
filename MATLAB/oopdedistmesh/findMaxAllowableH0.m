function max_h0 = findMaxAllowableH0(initial_h0, fd_circle, fh_non_uniform, bbox, pfix, geps, max_iterations, f, g, uD, N0, theta)
    % Initial binary search range
    h0_lower = 0.01; % Lower bound of h0
    h0_upper = 10 * initial_h0;  % Upper bound of h0
    tol = 1e-4;      % Tolerance for convergence

    % Define the distortion criteria threshold
    threshold_aspect_ratio = 6; % Adjust as needed

    % Binary search loop
    while h0_upper - h0_lower > tol
        % Choose midpoint
        h0 = (h0_lower + h0_upper) / 2;

        % Create initial mesh using distmesh
        [p_initial, ~, coordinates, elements] = do_distmesh2d_polyhedral(fd_circle, fh_non_uniform, h0, bbox, pfix, geps, max_iterations, []);

        % Call adaptiveMeshRef with current h0
        [refined_coordinates, refined_elements] = adaptiveMeshRef(N0, coordinates, elements, [], [], uD, computeBoundaries(elements), f, g);

        % Check if refined mesh meets distortion criteria
        meets_criteria = checkDistortion(refined_coordinates, refined_elements, threshold_aspect_ratio);

        % Update binary search range
        if meets_criteria
            % If meets criteria, increase h0
            h0_lower = h0;
        else
            % If doesn't meet criteria, decrease h0
            h0_upper = h0;
        end
    end

    % Final maximum allowable h0
    max_h0 = h0_lower;
end

% distortion criteria function
function meets_criteria = checkDistortion(coordinates, elements, threshold_aspect_ratio)
    % Compute aspect ratio for each element
    aspect_ratios = computeAspectRatios(coordinates, elements);

    % Check if all elements meet the aspect ratio criterion
    meets_criteria = all(aspect_ratios <= threshold_aspect_ratio);
end

function aspect_ratios = computeAspectRatios(coordinates, elements)
    % Compute aspect ratio for each element
    n_elements = size(elements, 1);
    aspect_ratios = zeros(n_elements, 1);

    for i = 1:n_elements
        % Get coordinates of vertices of the current element
        vertices = coordinates(elements(i, :), :);

        % Compute edge lengths of the triangle
        edge_lengths = [
            norm(vertices(2, :) - vertices(1, :)); % Edge 1-2
            norm(vertices(3, :) - vertices(2, :)); % Edge 2-3
            norm(vertices(1, :) - vertices(3, :))  % Edge 3-1
        ];

        % Compute aspect ratio of the triangle
        max_edge_length = max(edge_lengths);
        min_edge_length = min(edge_lengths);
        aspect_ratios(i) = max_edge_length / min_edge_length;
    end
end
