function adaptiveMeshRefinementstationary(coordinates, elements, dirichlet, neumann)
    % Parameters
    theta = 0.5; % Dörfler criterion parameter
    max_elements = 2811808; % Maximum number of elements
    tau = 0.03; % Tolerance for error estimator
    sigma = 0.25; % Coarsening parameter
    theta_refine = 0.25; % Refinement parameter

    % Initialize mesh and solution
    x = solveLaplace(coordinates, elements, dirichlet, neumann);

    % Initialize variables
    num_elements = size(elements, 1);
    mesh_history = cell(1, num_elements);
    error_history = zeros(1, num_elements);
    current_mesh = struct('coordinates', coordinates, 'elements', elements);

    % Adaptive algorithm loop
    while num_elements < max_elements
        % Compute error indicators
        eta_T = computeEtaR(x, current_mesh.coordinates, current_mesh.elements, dirichlet, neumann);

        % Check Dörfler criterion for refinement
        marked_elements = findMarkedElements(eta_T, theta);

        % Refine marked elements
        refined_mesh = refineMesh(current_mesh, marked_elements);

        % Compute solution on refined mesh
        x = solveLaplace(refined_mesh.coordinates, refined_mesh.elements, dirichlet, neumann);

        % Update variables
        num_elements = size(refined_mesh.elements, 1);
        mesh_history{num_elements} = refined_mesh;
        error_history(num_elements) = computeError(x, refined_mesh, tau);

        % Update current mesh
        current_mesh = refined_mesh;
    end

    % Plot results
    plotResults(mesh_history, error_history);
end

function x = solveLaplace(coordinates, elements, dirichlet, neumann)
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

function etaR = computeEtaR(x,coordinates,elements,dirichlet,neumann,f,g)
    [edge2nodes,element2edges,dirichlet2edges,neumann2edges] = provideGeometricData(elements,dirichlet,neumann);

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
    etaR = accumarray(element2edges(:), [dudn21; dudn32; dudn13], [size(edge2nodes,1) 1]);

    % Incorporate Neumann data
    if ~isempty(neumann)
        cn1 = coordinates(neumann(:,1),:);
        cn2 = coordinates(neumann(:,2),:);
        gmE = feval(g, (cn1 + cn2) / 2);
        etaR(neumann2edges) = etaR(neumann2edges) - sqrt(sum((cn2 - cn1).^2, 2)) .* gmE;
    end

    % Incorporate Dirichlet data
    etaR(dirichlet2edges) = 0;

    % Assemble edge contributions of indicators
    etaR = sum(etaR(element2edges).^2, 2);

    % Add volume residual to indicators
    fsT = feval(f, (c1 + (d21 + d31) / 3));
    etaR = etaR + (0.5 * area2 .* fsT).^2;
end

function marked_elements = findMarkedElements(eta_T, theta)
    [~, idx] = sort(eta_T, 'descend');
    sumeta = cumsum(eta_T);
    ell = find(sumeta >= sumeta(end) * theta, 1);
    marked_elements = idx(1:ell);
end

function refined_mesh = refineMesh(current_mesh, marked_elements)
    coordinates = current_mesh.coordinates;
    elements = current_mesh.elements;

    % Obtain geometric information on edges
    [edge2nodes, element2edges] = provideGeometricData(elements);

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

    % Refine boundary conditions
    boundary = cell(1, nargin - 2);
    for j = 1:numel(boundary)
        boundary{j} = [];
    end

    % Provide new nodes for refinement of elements
    newNodes = edge2newNode(element2edges);

    % Determine type of refinement for each element
    markedEdges = (newNodes ~= 0);
    none = ~markedEdges(:, 1);
    bisec1 = (markedEdges(:, 1) & ~markedEdges(:, 2) & ~markedEdges(:, 3));
    bisec12 = (markedEdges(:, 1) & markedEdges(:, 2) & ~markedEdges(:, 3));
    bisec13 = (markedEdges(:, 1) & ~markedEdges(:, 2) & markedEdges(:, 3));
    bisec123 = (markedEdges(:, 1) & markedEdges(:, 2) & markedEdges(:, 3));

    % Generate element numbering for refined mesh
    idx = ones(size(elements, 1), 1);
    idx(bisec1) = 2; % bisec(1): newest vertex bisection of 1st edge
    idx(bisec12) = 3; % bisec(2): newest vertex bisection of 1st and 2nd edge
    idx(bisec13) = 3; % bisec(2): newest vertex bisection of 1st and 3rd edge
    idx(bisec123) = 4; % bisec(3): newest vertex bisection of all edges
    idx = [1; 1 + cumsum(idx)];

    % Generate new elements
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

    refined_mesh = struct('coordinates', coordinates, 'elements', newElements);
end

function error = computeError(x, mesh, tau)
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
end

function plotResults(mesh_history, error_history)
    % Implement function to plot results
    % You can customize this function to visualize the results as needed.
end



