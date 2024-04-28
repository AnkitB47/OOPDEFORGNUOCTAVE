function [p, t] = do_distmesh2d_polyhedral(fd, fh, h0, bbox, pfix, geps, max_iterations, Fscale, varargin)
    dptol = 0.01;
    ttol = 0.1;
    Fscale = 1.2;
    deltat = 0.2;
    geps = 0.001 * h0;
    deps = sqrt(eps) * h0;
    densityctrlfreq = 30;
    maxIterations = 100;
    min_points_threshold = 5;
    N = 10;  % You can adjust this value based on your needs

    % Check if the vertices are 2D or 3D
    dimension = size(vertices, 2);

    % 1. Create initial distribution based on polyhedral geometry
    p = vertices(:, 1:dimension);

    % Initialize r0 before the loop
    r0 = zeros(size(p, 1), 1);

    remove_indices = find(feval(fd, p, varargin{:}) >= geps);
    fprintf('Size of p before removal: %d x %d\n', size(p, 1), size(p, 2));
    fprintf('Indices to be removed: %s\n', mat2str(remove_indices));
    p(remove_indices, :) = [];
    fprintf('Size of p after removal: %d x %d\n', size(p, 1), size(p, 2));



    % Check if p is empty before applying the rejection method
    if ~isempty(p)
        fh_output = feval(fh, p, varargin{:});
        disp(['Size of fh_output: ', mat2str(size(fh_output))]);
        % Ensure that fh_output is a column vector
        if size(fh_output, 2) > 1
            fh_output = fh_output';
        end
        % Display the size of r0 before the division
        disp(['Size of r0 before division: ', mat2str(size(r0))]);

        r0 = 1. / sum(fh_output.^2, 2);   % Probability to keep point
        r0 = r0(:);  % Ensure r0 is a column vector


        % Display the size of r0 after the division
        disp(['Size of r0 after division: ', mat2str(size(r0))]);
        fprintf('Size of r0: %d x %d\n', size(r0, 1), size(r0, 2));  % Debugging line
        fprintf('Max of r0: %f\n', max(r0));  % Debugging line
        fprintf('Size of rand(size(p, 1), 1): %d x %d\n', size(rand(size(p, 1), 1)));  % Debugging line
        fprintf('Size of max(r0): %d x %d\n', size(max(r0)));  % Debugging line

        % Perform rejection method
        max_rejection_attempts = 5;  % Set a maximum number of attempts
        rejection_attempt = 0;

        % Save the current points for debugging
        p_before_rejection = p;

       while rejection_attempt < max_rejection_attempts
       % Check if p is empty before applying the rejection method
            if ~isempty(p)
            % Evaluate the signed distance function at the current points
            signed_dist_before_rejection = feval(fd, p, varargin{:});

            % Display debugging information
            disp('Rejection method:');
            disp(['Before removal - Size of p: ', mat2str(size(p))]);
            disp(['Signed distances: ', mat2str(signed_dist_before_rejection')]);

            % Find indices of points violating the constraint
            indices_to_remove = find(signed_dist_before_rejection >= geps);

            % Display additional debugging information
            disp(['Indices to remove: ', mat2str(indices_to_remove')]);

            % Remove violating points
            if ~isempty(indices_to_remove)
                p(indices_to_remove, :) = [];

                % Update r0 after removal
                if ~isempty(p)
                    fh_output = feval(fh, p, varargin{:});
                    r0 = 1. / sum(fh_output.^2, 2);  % Probability to keep point
                else
                    r0 = [];
                end

                % Display some debugging information after removal
                disp(['After removal - Size of p: ', mat2str(size(p))]);

                % Break out of the while loop if removal was successful
                break;
            else
                disp('Error: All points were rejected. Retrying...');
            end
        else
            disp('Error: No points left to apply rejection method. Exiting loop.');
            break;
        end

        rejection_attempt = rejection_attempt + 1;

        if rejection_attempt >= max_rejection_attempts
            disp('Error: Maximum number of rejection attempts reached. Exiting loop.');
            break;
        end
    end

    % Display indices removed after the loop
    disp(['Indices removed: ', mat2str(indices_to_remove')]);

    end

    % Additional code for checking 2D or 3D
    if dimension == 2
        % Handle 2D case
        % Clean up the input points
        p = unique(p, 'rows', 'stable');  % Remove duplicate points
        p = p(all(isfinite(p), 2), :);     % Remove rows with NaN values

        try
            tri = delaunay(p(:, 1), p(:, 2));
            t = tri;
            e = [tri(:, [1, 2]); tri(:, [1, 3]); tri(:, [2, 3])];
        catch delaunay_error
            disp("Error in Delaunay triangulation in 2D:");
            disp(getReport(delaunay_error));
            t = [];
            e = [];
        end
    elseif dimension == 3
        % Handle 3D case
        % Clean up the input points
        p = unique(p, 'rows', 'stable');  % Remove duplicate points
        p = p(all(isfinite(p), 2), :);     % Remove rows with NaN values

        % Tetrahedral mesh generation in 3D (Delaunay)
        try
            tet = delaunayn(p);
            t = tet;
            e = [tet(:, [1, 2]); tet(:, [1, 3]); tet(:, [1, 4]); tet(:, [2, 3]); tet(:, [2, 4]); tet(:, [3, 4])];
        catch delaunayn_error
            disp("Error in Delaunay triangulation in 3D:");
            disp(getReport(delaunayn_error));
            t = [];
            e = [];
        end
    else
        error('Invalid dimension in the polyhedral geometry.');
    end

    % Initialize count outside of the loop
    count = 1;

    % Initialize pold before the loop
    pold = p;

    % Initialize bars outside of the loop
    bars = [];

    % Add a flag to indicate when to exit the loop
    exitLoop = false;

    node_count_history = zeros(maxIterations, 1);

    try
        while true
            disp(['Iteration: ', num2str(count), ' of ', num2str(maxIterations)]);
            disp(['Number of points: ', num2str(size(p, 1))]);
            disp(['Number of triangles: ', num2str(size(t, 1))]);
            disp(['Number of bars: ', num2str(size(bars, 1))]);

            if count > maxIterations
                disp('Maximum iterations reached. Exiting loop.');
                break;
            end

            % Add these lines for debugging
            disp('Debugging information:');
            disp(['pold size: ', mat2str(size(pold))]);
            disp(['p size: ', mat2str(size(p))]);
            disp(['t size: ', mat2str(size(t))]);
            disp(['bars size: ', mat2str(size(bars))]);

            % 3. Retriangulation by the Delaunay algorithm
            if count == 1 || max(sqrt(sum((p - pold).^2, 2)) / h0) > ttol  % Any large movement or first iteration?

                pold = p;  % Save current positions

                % Clean up the input points
                p = unique(p, 'rows', 'stable');  % Remove duplicate points
                p = p(all(isfinite(p), 2), :);     % Remove rows with NaN values

                try
                    % Using alternative triangulation code
                    t = delaunay(p(:, 1), p(:, 2));

                    % Check for duplicate triangles and remove them
                    [~, idx] = unique(sort(t, 2), 'rows');
                    t = t(idx, :);

                    % Remove triangles outside the region
                    pmid = (p(t(:, 1), :) + p(t(:, 2), :) + p(t(:, 3), :)) / 3;  % Compute centroids
                    interior_triangles = t(fd(pmid) < -geps, :);  % Keep interior triangles
                    t = interior_triangles;

                    % Print some debugging information
                    fprintf('Number of triangles before removal: %d\n', size(t, 1));
                    fprintf('Number of triangles after removal: %d\n', size(interior_triangles, 1));

                    % Print the indices of removed triangles
                    removed_indices = setdiff(1:size(t, 1), find(fd(pmid) < -geps));
                    fprintf('Indices of removed triangles: %s\n', mat2str(removed_indices));

                    % Describe each bar by a unique pair of nodes
                    if ~isempty(t)
                        edges = [t(:, [1, 2]); t(:, [1, 3]); t(:, [2, 3])];  % Interior bars duplicated
                        edges = unique(sort(edges, 2), 'rows');  % Remove duplicate edges

                        % Ensure that edges does not exceed N rows
                        edges = edges(1:min(end, N), :);

                        % Update the bars list
                        bars = edges;

                        % Print some debugging information about bars
                        fprintf('Number of bars before removal: %d\n', size(bars, 1));

                        % Remove duplicate bars
                        [~, idx] = unique(sort(bars, 2), 'rows');
                        bars = bars(idx, :);

                        % Ensure that bars does not exceed N rows
                        bars = bars(1:min(end, N), :);

                        % Print some debugging information about bars after removal
                        fprintf('Number of bars after removal: %d\n', size(bars, 1))

                    else
                        disp('Error: No valid triangles in the mesh.');
                        exitLoop = true;  % Set flag to exit the loop
                    end

                catch triangulation_error
                    disp("Error in triangulation in 2D:");
                    disp(getfield(triangulation_error, 'message'));
                    t = [];
                end
            end

            % Check if bars is not empty before proceeding
            if ~isempty(bars)
                % 6. Move mesh points based on bar lengths L and forces F
                barvec = p(bars(:, 1), :) - p(bars(:, 2), :);
                L = sqrt(sum(barvec.^2, 2));
                hbars = fh((p(bars(:, 1), :) + p(bars(:, 2), :)) / 2);
                L0 = hbars * Fscale * sqrt(sum(L.^2) / sum(hbars.^2));  % L0 = Desired lengths
                F = (1 - L0 ./ L) / 2;  % Nodal forces

                Fj1 = sparse(bars(:, 1), ones(size(bars, 1), 1), F, size(p, 1), 1);
                Fj2 = sparse(bars(:, 2), ones(size(bars, 1), 1), F, size(p, 1), 1);
                F = Fj1 + Fj2;  % Combine forces

                % Updated code for better readability and debugging
                % Compute the magnitudes of forces
                magnitude_F = sqrt(sum(F.^2, 2));

                % Ensure that the magnitudes are not zero to avoid division by zero
                magnitude_F(magnitude_F == 0) = 1;  % Set to 1 to avoid division by zero

                % Normalize the forces
                normalized_forces = F ./ [magnitude_F, magnitude_F];

                % Accumulate forces for each node
                force_accumulator = zeros(size(p));
                for k = 1:size(bars, 1)
                    force_accumulator(bars(k, 1), :) = force_accumulator(bars(k, 1), :) + normalized_forces(k, :);
                    force_accumulator(bars(k, 2), :) = force_accumulator(bars(k, 2), :) + normalized_forces(k, :);
                end

                % Update node positions
                p = p + deltat * force_accumulator .* p;

                % 7. Move mesh points based on bar forces F
                p = p + deltat * F ./ (sqrt(sum(F.^2, 2)) * [1, 1]) .* p;  % Update node positions

                if ~isempty(pfix)
                    p(1:nfix, :) = pfix;
                end  % Reset fixed points

                % 8. Bring outside points back to the boundary
                p(p(:, 1) < bbox(1, 1), 1) = bbox(1, 1); p(p(:, 1) > bbox(2, 1), 1) = bbox(2, 1);
                p(p(:, 2) < bbox(1, 2), 2) = bbox(1, 2); p(p(:, 2) > bbox(2, 2), 2) = bbox(2, 2);

                % 9. Terminate if the mesh has converged
                node_count = size(p, 1);
                node_count_history(count) = node_count;

                % Check for convergence based on the change in the number of nodes
                if count > 5
                    % Calculate the relative change in node count over the last 5 iterations
                    rel_change = abs(diff(node_count_history(count-4:count-1))) / node_count;

                    % If the maximum relative change is less than a threshold, consider the mesh converged
                    if max(rel_change) < 1e-3
                        disp('Mesh has converged. Exiting loop.');
                        break;
                    end
                end

                % 10. Mesh density control
                if rem(count, densityctrlfreq) == 0
                    % Density control code
                    % Calculate local density at each node
                    node_density = zeros(size(p, 1), 1);
                    for k = 1:size(bars, 1)
                        node_density(bars(k, 1)) = node_density(bars(k, 1)) + 1;
                        node_density(bars(k, 2)) = node_density(bars(k, 2)) + 1;
                    end

                    % Evaluate mesh density function fdensity(p) at each node
                    density_values = feval(fdensity, p);

                    % Update node positions based on density control
                    p = p + deltat_density * density_values .* (pfix - p);
                end
                count = count + 1;
            end

            % Exit the loop if the flag is set
            if exitLoop
                break;
            end
        end
        disp('Termination conditions met. Exiting loop.');
        disp('Mesh generation completed.');
    catch exception
        disp('Error occurred:');
        disp(exception.message);
    end
end

function force_accumulator = accumulateForces(bars, F, N)
    force_accumulator = zeros(N, 2);
    for k = 1:size(bars, 1)
        force_accumulator(bars(k, 1), :) = force_accumulator(bars(k, 1), :) + F(k, :);
        force_accumulator(bars(k, 2), :) = force_accumulator(bars(k, 2), :) + F(k, :);
    end
end

function p = enforceBoundary(p, bbox)
    p(p(:, 1) < bbox(1, 1), 1) = bbox(1, 1);
    p(p(:, 1) > bbox(2, 1), 1) = bbox(2, 1);
    p(p(:, 2) < bbox(1, 2), 2) = bbox(1, 2);
    p(p(:, 2) > bbox(2, 2), 2) = bbox(2, 2);
end

function node_density = accumulateDensity(bars, N)
    node_density = zeros(N, 1);
    for k = 1:size(bars, 1)
        node_density(bars(k, 1)) = node_density(bars(k, 1)) + 1;
        node_density(bars(k, 2)) = node_density(bars(k, 2)) + 1;
    end
end

function p = densityControl(density_values, pfix, p)
    p = p + density_values .* (pfix - p);
end
