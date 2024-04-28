function [coordinates, elements] = adaptiveMeshRef(N0, coordinates, elements, dirichlet, neumann, uD, boundaries, f, g)
    % Copyright (C) 2024 Ankit Bhardwaj
    % Parameters
    theta = 0.5; % Dörfler criterion parameter
    max_elements = 10000; % Maximum number of elements
    tau = 1e-6; % Tolerance for error estimator
    time_max = 1e+6;

    % Initialize variables
    num_elements = size(elements, 1);
    solve_brute_cumulative_time = zeros(1, time_max);
    solve_optimised_cumulative_time = zeros(1, time_max);
    refine_cumulative_time = zeros(1, time_max);
    coarsen_cumulative_time = zeros(1, time_max);
    solve_brute_time = 0;
    solve_optimised_time = 0;
    refine_time = 0;
    coarsen_time = 0;

    % Adaptive algorithm loop
    while num_elements < max_elements
        % Compute solution time for brute version
        tic;
        x_brute = solveLaplace_brute(coordinates, elements, dirichlet, neumann, uD, f, g);
        solve_brute_time = solve_brute_time + toc;
        solve_brute_cumulative_time(num_elements) = solve_brute_time;

        % Compute solution time for optimised version
        tic;
        x_optimised = solveLaplace_optimised(coordinates, elements, dirichlet, neumann, uD, f, g);
        solve_optimised_time = solve_optimised_time + toc;
        solve_optimised_cumulative_time(num_elements) = solve_optimised_time;

        % Compute error indicators
        eta_T = computeEtaR(x_optimised, coordinates, elements, dirichlet, neumann, boundaries, f, g);

        % Check Dörfler criterion for refinement
        marked_elements = findMarkedElements(eta_T, theta);

        % Refine mesh and measure time
        tic;
        [coordinates, elements] = refineMesh(coordinates, elements, marked_elements, boundaries);
        refine_time = refine_time + toc;
        refine_cumulative_time(num_elements) = refine_time;

        % Coarsen mesh and measure time
        tic;
        [coordinates, elements] = coarsenNVB(N0, coordinates, elements, marked_elements, boundaries);
        coarsen_time = coarsen_time + toc;
        coarsen_cumulative_time(num_elements) = coarsen_time;

        % Update variables
        num_elements = size(elements, 1);

        % Check if the maximum residual error is less than the tolerance
        if max(eta_T) < tau
            disp('Tolerance reached. Terminating...');
            break;
        end

        if num_elements >= max_elements
            disp('Number of elements close to the maximum allowed. Terminating...');
            break;
        end
    end

    % Plot results
    figure;

    % Plot brute solveLaplace
    subplot(1, 2, 1);
    plot(1:num_elements, solve_brute_cumulative_time(1:num_elements)*1000, '-o', 'DisplayName', 'solveLaplace Brute');
    ylabel('Cumulative Computational Time (milliseconds)');
    xlabel('Number of Elements (M)');
    title('Brute solveLaplace');
    grid on;
    xlim([1, 10000]);

    % Plot other methods
    subplot(1, 2, 2);
    plot(1:num_elements, solve_optimised_cumulative_time(1:num_elements)*1000, '-o', 'DisplayName', 'solveLaplace Optimised');
    hold on;
    plot(1:num_elements, refine_cumulative_time(1:num_elements)*1000, '-o', 'DisplayName', 'refineMesh');
    plot(1:num_elements, coarsen_cumulative_time(1:num_elements)*1000, '-o', 'DisplayName', 'coarsenNVB');
    hold off;
    ylabel('Cumulative Computational Time (milliseconds)');
    xlabel('Number of Elements (M)');
    title('Other Methods');
    grid on;
    xlim([1, 10000]);

    legend;

end

function [edge2nodes, element2edges, dirichlet2edges, neumann2edges] = provideGeometricData(elements, boundaries)
    nE = size(elements, 1);
    nB = length(boundaries);

    % Node vectors of all edges (interior edges appear twice)
    I = elements(:);
    J = reshape(elements(:, [2, 3, 1]), 3 * nE, 1);

    % Create edge matrix
    edges = sort([I J], 2);

    % Find unique edges and their indices
    [uniqueEdges, ~, edgeNumber] = unique(edges, 'rows');

    % Provide element2edges and edge2nodes
    element2edges = reshape(edgeNumber, nE, 3);
    edge2nodes = uniqueEdges;

    % Provide boundary2edges
    boundary2edges = cell(1, nB);
    for j = 1:nB
        boundary = boundaries{j};
        if ~isempty(boundary)
            boundaryEdges = sort(boundary, 2);
            boundaryIndices = zeros(size(boundaryEdges, 1), 1);
            for k = 1:size(boundaryEdges, 1)
                idx = find(ismember(uniqueEdges, boundaryEdges(k, :), 'rows'));
               if ~isempty(idx)
                    boundaryIndices(k) = idx;
                end
            end
            boundary2edges{j} = boundaryIndices;
        else
            boundary2edges{j} = [];
        end
    end

    % Print the dimensions of element2edges and edge2nodes
    fprintf('Size of element2edges: %dx%d\n', size(element2edges, 1), size(element2edges, 2));
    fprintf('Size of edge2nodes: %dx%d\n', size(edge2nodes, 1), size(edge2nodes, 2));

    % Print a few values of element2edges and edge2nodes for inspection
    disp('element2edges:');
    disp(element2edges(1:min(5, size(element2edges, 1)), :));
    disp('edge2nodes:');
    disp(edge2nodes(1:min(5, size(edge2nodes, 1)), :));

    % Dummy output for demonstration
    dirichlet2edges = [];
    neumann2edges = [];
end

function A = assembleStiffnessMatrix_LinearComplexity(coordinates, elements)
    nE = size(elements, 1);
    nC = size(coordinates, 1);

    I = zeros(9 * nE, 1);
    J = zeros(9 * nE, 1);
    A = zeros(9 * nE, 1);

    for i = 1:nE
        nodes = elements(i,:);
        B = [1, 1, 1; coordinates(nodes,:)'];
        grad = B \ [0, 0; 1, 0; 0, 1];
        idx = 9 * (i - 1) + 1:9 * i;
        tmp = [1; 1; 1] * nodes;
        I(idx) = reshape(tmp', 9, 1);
        J(idx) = reshape(tmp, 9, 1);
        A(idx) = det(B) / 2 * reshape(grad * grad', 9, 1);
    end

    A = sparse(I, J, A, nC, nC);
end

function [x, energy] = solveLaplace_brute(coordinates, elements, dirichlet, neumann, f, g, uD)
    nC = size(coordinates, 1);
    x = zeros(nC, 1);

    % Assembly of stiffness matrix
    A = sparse(nC, nC);
    for i = 1:size(elements, 1)
        nodes = elements(i,:);
        B = [1, 1, 1; coordinates(nodes,:)'];
        grad = B \ [0, 0; 1, 0; 0, 1];
        A(nodes, nodes) = A(nodes, nodes) + det(B) * grad * grad' / 2;
    end

    % Prescribe values at Dirichlet nodes
    dirichlet = unique(dirichlet);
    x(dirichlet) = feval(uD, coordinates(dirichlet,:));

    % Assembly of right-hand side
    b = -A * x;
    for i = 1:size(elements, 1)
        nodes = elements(i,:);
        sT = [1, 1, 1] * coordinates(nodes,:) / 3;
        b(nodes) = b(nodes) + det([1, 1, 1; coordinates(nodes,:)']) * feval(f, sT) / 6;
    end
    for i = 1:size(neumann, 1)
        nodes = neumann(i,:);
        mE = [1, 1] * coordinates(nodes,:) / 2;
        b(nodes) = b(nodes) + norm([1, -1] * coordinates(nodes,:)) * feval(g, mE) / 2;
    end

    % Computation of P1-FEM approximation
    freenodes = setdiff(1:nC, dirichlet);
    x(freenodes) = A(freenodes, freenodes) \ b(freenodes);

    % Compute energy | | grad(uh) | | ˆ 2 of discrete solution
    energy = x' * A * x;
end

function x = solveLaplace_optimised(coordinates, elements, dirichlet, neumann, uD, f, g)
    nE = size(elements, 1);
    nC = size(coordinates, 1);
    x = zeros(nC, 1);

    % First vertex of elements and corresponding edge vectors
    c1 = coordinates(elements(:, 1), :);
    d21 = coordinates(elements(:, 2), :) - c1;
    d31 = coordinates(elements(:, 3), :) - c1;

    % Vector of element areas 4*|T|
    area4 = 2 * (d21(:, 1) .* d31(:, 2) - d21(:, 2) .* d31(:, 1));

    % Assembly of stiffness matrix
    I = reshape(elements(:, [1 2 3 1 2 3 1 2 3])', 9 * nE, 1);
    J = reshape(elements(:, [1 1 1 2 2 2 3 3 3])', 9 * nE, 1);
    a = (sum(d21 .* d31, 2) ./ area4)';
    b = (sum(d31 .* d31, 2) ./ area4)';
    c = (sum(d21 .* d21, 2) ./ area4)';
    A = [-2 * a + b + c; a - b; a - c; a - b; b; -a; a - c; -a; c];
    A = sparse(I, J, A(:));

    % Prescribe values at Dirichlet nodes
    dirichlet = unique(dirichlet);
    x(dirichlet) = uD(coordinates(dirichlet, :));

    % Assembly of right-hand side
    fsT = f(c1 + (d21 + d31) / 3);
    b = accumarray(elements(:), repmat(12 \ area4 .* fsT, 3, 1), [nC 1]) - A * x;
    if ~isempty(neumann)
        cn1 = coordinates(neumann(:, 1), :);
        cn2 = coordinates(neumann(:, 2), :);
        gmE = g((cn1 + cn2) / 2);
        b = b + accumarray(neumann(:), repmat(2 \ sqrt(sum((cn2 - cn1).^2, 2)) .* gmE, 2, 1), [nC 1]);
    end

    % Computation of P1-FEM approximation
    freenodes = setdiff(1:nC, dirichlet);
    x(freenodes) = A(freenodes, freenodes) \ b(freenodes);
    energy = x' * A * x;
end

function marked_elements = findMarkedElements(eta_T, theta)
    [~, idx] = sort(eta_T, 'descend');
    sumeta = cumsum(eta_T);
    ell = find(sumeta >= sumeta(end) * theta, 1);
    marked_elements = idx(1:ell);
end

function etaR = computeEtaR(x, coordinates, elements, dirichlet, neumann, boundaries, f, g)
    [edge2nodes, element2edges, dirichlet2edges, neumann2edges] = provideGeometricData(elements, boundaries);

    % First vertex of elements and corresponding edge vectors
    c1 = coordinates(elements(:,1),:);
    d21 = coordinates(elements(:,2),:) - c1;
    d31 = coordinates(elements(:,3),:) - c1;

    % Vector of element volumes 2*|T |
    area2 = d21(:,1) .* d31(:,2) - d21(:,2) .* d31(:,1);

    % Compute curl(uh) = (-duh/dy, duh/dx)
    u21 = repmat(x(elements(:,2)) - x(elements(:,1)), 1, 2);
    u31 = repmat(x(elements(:,3)) - x(elements(:,1)), 1, 2);
    curl = (d31 .* u21 - d21 .* u31) ./ repmat(area2, 1, 2);

    % Compute edge terms hE*(duh/dn) for uh
    dudn21 = sum(d21 .* curl, 2);
    dudn13 = -sum(d31 .* curl, 2);
    dudn32 = -(dudn13 + dudn21);

    % Compute residual-based error estimator
    etaR = sqrt(dudn21.^2 + dudn32.^2 + dudn13.^2);

    % Incorporate Neumann data
    if ~isempty(neumann)
        cn1 = coordinates(neumann(:,1),:);
        cn2 = coordinates(neumann(:,2),:);
        gmE = feval(g, (cn1 + cn2) / 2); % Compute g at midpoints
        lengths = sqrt(sum((cn2 - cn1).^2, 2)); % Compute lengths of Neumann edges
        etaR(neumann2edges) = etaR(neumann2edges) - lengths .* gmE;
    end

    % Incorporate Dirichlet data
    etaR(dirichlet2edges) = 0;

    % Assemble edge contributions of indicators
    etaR = sum(etaR.^2, 2);

    % Add volume residual to indicators
    fsT = feval(f, (c1 + (d21 + d31) / 3)); % Compute f at centroids
    etaR = etaR + (0.5 * area2 .* fsT).^2;

    % Compute residual errors
    residualErrors = 0.5 * area2 .* fsT;

    % Display residual errors
    disp('Residual Errors:');
    disp(residualErrors);

    % Display residual-based error estimator
    disp('Residual-based Error Estimator:');
    disp(etaR);

end

function [coordinates, elements] = refineMesh(coordinates, elements, marked_elements, boundaries)

    disp('Inside refineMesh function.');

    % Integrate provideGeometricData function
    [edge2nodes, element2edges,~, ~] = provideGeometricData(elements, boundaries);

    % Mark edges for refinement
    edge2newNode = zeros(max(max(element2edges)), 1);
    edge2newNode(element2edges(marked_elements, :)) = 1;
    swap = 1;
    while ~isempty(swap)
        markedEdge = edge2newNode(element2edges);
        swap = find(~markedEdge(:, 1) & (markedEdge(:, 2) | markedEdge(:, 3)));
        edge2newNode(element2edges(swap, 1)) = 1;
    end

    % Generate new nodes
    edge2newNode(edge2newNode ~= 0) = size(coordinates, 1) + (1:nnz(edge2newNode));
    idx = find(edge2newNode);
    coordinates(edge2newNode(idx), :) = (coordinates(edge2nodes(idx, 1), :) + coordinates(edge2nodes(idx, 2), :)) / 2;

    % Refine elements
    newNodes = edge2newNode(element2edges);
    markedEdges = (newNodes ~= 0);
    none = ~markedEdges(:, 1);
    bisec1 = (markedEdges(:, 1) & ~markedEdges(:, 2) & ~markedEdges(:, 3));
    bisec12 = (markedEdges(:, 1) & markedEdges(:, 2) & ~markedEdges(:, 3));
    bisec13 = (markedEdges(:, 1) & ~markedEdges(:, 2) & markedEdges(:, 3));
    bisec123 = (markedEdges(:, 1) & markedEdges(:, 2) & markedEdges(:, 3));
    idx = ones(size(elements, 1), 1);
    idx(bisec1) = 2;
    idx(bisec12) = 3;
    idx(bisec13) = 3;
    idx(bisec123) = 4;
    idx = [1; 1 + cumsum(idx)];
    newElements = zeros(idx(end) - 1, 3);
    newElements(idx(none), :) = elements(none, :);
    newElements([idx(bisec1), 1 + idx(bisec1)], :) = [elements(bisec1, 3), elements(bisec1, 1), newNodes(bisec1, 1); ...
                                                     elements(bisec1, 2), elements(bisec1, 3), newNodes(bisec1, 1)];
    newElements([idx(bisec12), 1 + idx(bisec12), 2 + idx(bisec12)], :) = [elements(bisec12, 3), elements(bisec12, 1), newNodes(bisec12, 1); ...
                                                                          newNodes(bisec12, 1), elements(bisec12, 2), newNodes(bisec12, 2); ...
                                                                          elements(bisec12, 3), newNodes(bisec12, 1), newNodes(bisec12, 2)];
    newElements([idx(bisec13), 1 + idx(bisec13), 2 + idx(bisec13)], :) = [newNodes(bisec13, 1), elements(bisec13, 3), newNodes(bisec13, 3); ...
                                                                          elements(bisec13, 1), newNodes(bisec13, 1), newNodes(bisec13, 3); ...
                                                                          elements(bisec13, 2), elements(bisec13, 3), newNodes(bisec13, 1)];
    newElements([idx(bisec123), 1 + idx(bisec123), 2 + idx(bisec123), 3 + idx(bisec123)], :) = [newNodes(bisec123, 1), elements(bisec123, 3), newNodes(bisec123, 3); ...
                                                                                                  elements(bisec123, 1), newNodes(bisec123, 1), newNodes(bisec123, 3); ...
                                                                                                  newNodes(bisec123, 1), elements(bisec123, 2), newNodes(bisec123, 2); ...
                                                                                                  elements(bisec123, 3), newNodes(bisec123, 1), newNodes(bisec123, 2)];

    % Update coordinates and elements
    coordinates = coordinates;
    elements = newElements;

    % Display the size of updated coordinates and elements
    disp(['Size of updated coordinates: ', num2str(size(coordinates))]);
    disp(['Size of updated elements: ', num2str(size(elements))]);

end

function [coordinates, elements, varargout] = coarsenNVB(N0, coordinates, elements, varargin)
    nC = size(coordinates, 1);
    nE = size(elements, 1);

    % Obtain geometric information on neighbouring elements
    I = elements(:);
    J = reshape(elements(:, [2, 3, 1]), 3 * nE, 1);
    nodes2edge = sparse(I, J, 1:3 * nE);
    mask = nodes2edge > 0;
    [foo{1:2}, idxIJ] = find(nodes2edge);
    [foo{1:2}, neighbourIJ] = find(mask + mask .* sparse(J, I, [1:nE, 1:nE, 1:nE]'));
    element2neighbours(idxIJ) = neighbourIJ - 1;
    element2neighbours = reshape(element2neighbours, nE, 3);

    % Determine which nodes are deleted by coarsening
    marked = zeros(nE, 1);
    marked(varargin{end}{1}) = 1;
    newestNode = unique(elements((marked & elements(:, 3) > N0), 3));
    valence = accumarray(elements(:), 1, [nC, 1]);
    markedNodes = zeros(nC, 1);
    markedNodes(newestNode((valence(newestNode) == 2 | valence(newestNode) == 4))) = 1;

    % Collect pairs of brother elements that will be united
    idx = find(markedNodes(elements(:, 3)) & (element2neighbours(:, 3) > (1:nE)'))';
    markedElements = zeros(nE, 1);
    markedElements(idx) = 1;
    for element = idx
        if markedElements(element)
            markedElements(element2neighbours(element, 3)) = 0;
        end
    end
    idx = find(markedElements);

    % Coarsen two brother elements
    brother = element2neighbours(idx, 3);
    elements(idx, [1 3 2]) = [elements(idx, [2 1]) elements(brother, 1)];

    % Delete redundant nodes
    activeNodes = find(~markedNodes);
    coordinates = coordinates(activeNodes, :);

    % Provide permutation of nodes to correct further data
    coordinates2newCoordinates = zeros(1, nC);
    coordinates2newCoordinates(activeNodes) = 1:length(activeNodes);

    % Delete redundant elements + correct elements
    elements(brother, :) = [];
    elements = coordinates2newCoordinates(elements);

    % Delete redundant boundaries + correct boundaries
    for j = 1:nargout - 2
        boundary = varargin{j};
        if ~isempty(boundary)
            node2boundary = zeros(nC, 2);
            node2boundary(boundary(:, 1), 1) = 1:size(boundary, 1);
            node2boundary(boundary(:, 2), 2) = 1:size(boundary, 1);
            idx = (markedNodes & node2boundary(:, 2));
            boundary(node2boundary(idx, 2), 2) = boundary(node2boundary(idx, 1), 2);
            boundary(node2boundary(idx, 1), 2) = 0;
            varargout{j} = coordinates2newCoordinates(boundary(find(boundary(:, 2)), :));
        else
            varargout{j} = [];
        end
    end
end

function [x, coordinates, elements, indicators] = adaptiveAlgorithm(coordinates, elements, dirichlet, neumann, f, g, uD, nEmax, theta)
    while 1
        x = solveLaplace(coordinates, elements, dirichlet, neumann, f, g, uD);
        indicators = computeEtaR(x, coordinates, elements, dirichlet, neumann, f, g);

        if size(elements,1) >= nEmax
            break
        end

        [indicators,idx] = sort(indicators,'descend');
        sumeta = cumsum(indicators);
        ell = find(sumeta >= sumeta(end) * theta, 1);
        marked = idx(1:ell);

        [coordinates, elements, dirichlet, neumann] = refineNVB(coordinates, elements, dirichlet, neumann, marked);
    end
end

