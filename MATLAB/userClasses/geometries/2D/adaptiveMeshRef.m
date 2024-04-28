function [coordinates, elements] = adaptiveMeshRef(N0, coordinates, elements, dirichlet, neumann, uD, boundaries, f, g)
    % Parameters
    theta = 0.5; % Dörfler criterion parameter
    max_elements = 1000; % Maximum number of elements
    tau = 0.03; % Tolerance for error estimator

    % Initialize mesh and solution
    x = solveLaplace(coordinates, elements, dirichlet, neumann, uD, f, g);

    % Initialize variables
    num_elements = size(elements, 1);
    mesh_history = cell(1, num_elements);
    error_history = zeros(1, num_elements);

    % Adaptive algorithm loop
    while num_elements < max_elements
        % Compute error indicators
        eta_T = computeEtaR(x, coordinates, elements, dirichlet, neumann, boundaries, f, g);

        % Check Dörfler criterion for refinement
        marked_elements = findMarkedElements(eta_T, theta);

        % Refine mesh
        [coordinates, elements] = refineMesh(coordinates, elements, marked_elements, boundaries);

        % Coarsen mesh
        [coordinates, elements] = coarsenNVB(N0, coordinates, elements, marked_elements, boundaries);

        % Compute solution on refined mesh
        x = solveLaplace(coordinates, elements, dirichlet, neumann, uD, f, g);

        % Update variables
        num_elements = size(elements, 1);
        mesh_history{num_elements} = struct('coordinates', coordinates, 'elements', elements);
        error_history(num_elements) = computeError(x, struct('coordinates', coordinates, 'elements', elements), tau);
    end

    % Plot results
    %plotResults(mesh_history, error_history);
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

function x = solveLaplace(coordinates, elements, dirichlet, neumann, uD, f, g)
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

    % Check if all arrays have the same number of elements
    if numel(dudn21) ~= numel(dudn32) || numel(dudn21) ~= numel(dudn13)
        error('Dimensions of filtered arrays do not match.');
    end

    % Filter out invalid indices
    valid_indices = (dudn21 ~= 0) | (dudn32 ~= 0) | (dudn13 ~= 0);
    dudn21 = dudn21(valid_indices);
    dudn32 = dudn32(valid_indices);
    dudn13 = dudn13(valid_indices);

    % Find unique rows in element2edges
    [~, ~, idx] = unique(element2edges, 'rows');

    % Compute valid_indices
    valid_indices = idx(1:numel(idx));

    % Compute indices using linear indexing
    max_index = size(edge2nodes, 1);
    num_elements = size(element2edges, 1);
    valid_element2edges = element2edges(valid_indices, :);
    valid_indices = (1:numel(valid_indices))';

    disp(size(valid_indices));

    % Compute etaR
    etaR = zeros(num_elements, 1);
    for i = 1:num_elements
        row = valid_element2edges(i, 1);
        col = valid_element2edges(i, 2);
        if row <= size(dudn21, 1) && col <= size(dudn21, 2)
            etaR(i) = sqrt(dudn21(row, col)^2 + dudn32(row, col)^2 + dudn13(row, col)^2);
        else
            etaR(i) = 0;
        end
    end

    % Incorporate Neumann data
    if ~isempty(neumann)
        cn1 = coordinates(neumann(:,1),:);
        cn2n2n2 = coordinates(neumann(:,2),:);
        gmE = feval(g, (cn1 + cn2)) = coordinates(neumann / 2);
        etaR(neumann2edges) = etaR(neumann2edges) - sqrt(sum((cn2 - cn1).^2, 2)) .* gmE;
    end

    coordinates(neumann) % Incorporate Dirichlet data
    etaR(dirichlet2edges) = 0;

    % Incorporate Dirichlet data
    etaR(dirichlet2edges) = 0;

    % Assemble edge contributions of indicators
    etaR = sum(etaR(valid_indices).^2, 2);

    disp(size(etaR));

    % Add volume residual to indicators
    fsT = feval(f, (c1 + (d21 + d31) / 3));
    etaR = etaR + (0.5 * area2 .* fsT).^2;

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

function error = computeError(x, mesh, tau)

    disp('Inside computeError function.');

    coordinates = mesh.coordinates;
    elements = mesh.elements;

    % First vertex of elements and corresponding edge vectors
    c1 = coordinates(elements(:,1),:);
    d21 = coordinates(elements(:,2),:) - c1;
    d31 = coordinates(elements(:,3),:) - c1;

    % Vector of element volumes 2*|T |
    area2 = d21(:,1) .* d31(:,2) - d21(:,2) .* d31(:,1);

    % Compute gradient of x
    grad_x = [-(d21(:, 2) .* x(elements(:,3)) - d31(:, 2) .* x(elements(:,2))) ./ area2, ...
              (d21(:, 1) .* x(elements(:,3)) - d31(:, 1) .* x(elements(:,2))) ./ area2];

    % Compute error
    error = sum(area2 .* (sum(grad_x.^2, 2) / 2)).^0.5;

    % Apply tolerance threshold
    error = max(error, tau);

    % Display the computed error
    disp(['Computed error: ', num2str(error)]);
end

function plotResults(mesh_history, error_history)
    % Check if mesh_history and error_history are empty
    if isempty(mesh_history) || isempty(error_history)
        error('mesh_history or error_history is empty. Cannot plot results.');
    end

    % Plot the history of mesh elements and corresponding error indicators
    figure;
    subplot(2, 1, 1);
    num_elements = zeros(size(mesh_history));
    for i = 1:numel(mesh_history)
        num_elements(i) = size(mesh_history{i}.elements, 1);
    end
    plot(1:numel(mesh_history), num_elements, 'b.-');
    xlabel('Iteration');
    ylabel('Number of Elements');
    title('Mesh Refinement History');

    subplot(2, 1, 2);
    plot(1:numel(error_history), error_history, 'r.-');
    xlabel('Iteration');
    ylabel('Error Indicator');
    title('Error Indicator History');
end


