function [x, coordinates, elements, indicators] = adaptiveMeshRefinement(coordinates, elements, dirichlet, neumann, f, g, uD, tmax, nEmax, theta, tau, sigma)
    n = 0;
    t = 0;
    delta_t = 0.1;  % Initialize delta_t with a value (adjust as needed)

    while t < tmax
        k = 0;
        Tn = elements;

        % Adaptive mesh refinement loop
        while true
            k = k + 1;

            % Compute discrete solution on current mesh
            x = solveLaplace(coordinates, Tn, dirichlet, neumann, f, g, uD);

            % Compute refinement indicators
            indicators = computeEtaR(x, coordinates, Tn, dirichlet, neumann, f, g);

            % Stopping criterion
            if size(Tn, 1) >= nEmax
                break;
            end

            % Mark elements for refinement
            [indicators, idx] = sort(indicators, 'descend');
            sumeta = cumsum(indicators);
            ell = find(sumeta >= sumeta(end) * theta, 1);
            marked = idx(1:ell);

            % Refine mesh using NVB
            [coordinates, Tn, dirichlet, neumann] = refineNVB(coordinates, Tn, dirichlet, neumann, marked);
        end

        % Use Dörfler criterion to mark elements for coarsening
        markedCoarse = findEtaCoarse(x, coordinates, Tn, dirichlet, neumann, f, g, sigma, tau);

        % Adaptive mesh coarsening loop
        while k > 1 && ~isempty(markedCoarse)
            k = k - 1;

            % Coarsen mesh
            [coordinates, Tn, dirichlet, neumann] = coarsenNVB(0, coordinates, Tn, dirichlet, neumann, markedCoarse);

            % Compute discrete solution on current mesh
            x = solveLaplace(coordinates, Tn, dirichlet, neumann, f, g, uD);

            % Compute refinement indicators
            indicators = computeEtaR(x, coordinates, Tn, dirichlet, neumann, f, g);

            % Mark elements for coarsening
            markedCoarse = findEtaCoarse(x, coordinates, Tn, dirichlet, neumann, f, g, sigma, tau);
        end

        % Go to next time step
        n = n + 1;
        t = t + delta_t;  % Assuming delta_t is defined
    end
end

function [edge2nodes, element2edges, dirichlet2edges, neumann2edges] = provideGeometricData(elements, dirichlet, neumann)
    % Function to provide geometric data for a triangular mesh

    % Get edges and nodes information
    edges = [elements(:, [1, 2]); elements(:, [2, 3]); elements(:, [3, 1])];
    nodes = unique(edges(:));

    % Create maps from edges to nodes
    edge2nodes = zeros(max(nodes), 3);
    for i = 1:size(edges, 1)
        edge2nodes(edges(i, :), mod(i - 1, 3) + 1) = 1;
    end

    % Create maps from elements to edges
    element2edges = zeros(size(elements, 1), max(nodes));
    for i = 1:size(elements, 1)
        element2edges(i, elements(i, :)) = 1;
    end

    % Ensure dirichlet values are within valid range
    validDirichlet = dirichlet(dirichlet >= 1 & dirichlet <= size(elements, 1));

    % Create dirichlet2edges with size equal to the number of nodes
    dirichlet2edges = zeros(max(nodes), max(nodes));

    % Check if dirichlet is not empty
    if ~isempty(validDirichlet)
        % Assign edges corresponding to Dirichlet boundary to dirichlet2edges
        for i = 1:numel(validDirichlet)
            elementIdx = validDirichlet(i);
            edgesIdx = find(element2edges(elementIdx, :));
            dirichlet2edges(edges(edgesIdx, 1), edges(edgesIdx, 2)) = 1;
            dirichlet2edges(edges(edgesIdx, 2), edges(edgesIdx, 1)) = 1;
        end
    end

    % Check if neumann is empty
    if ~isempty(neumann)
        % Create maps from Neumann edges to edges
        neumannEdges = unique(sort(neumann, 2), 'rows');
        neumann2edges = zeros(max(nodes), max(nodes));
        neumann2edges(neumannEdges(:, 1), neumannEdges(:, 2)) = 1;
    else
        neumann2edges = [];
    end
end

function [x, energy] = solveLaplace(coordinates, elements, dirichlet, neumann, f, g, uD)
    nE = size(elements, 1);
    nC = size(coordinates, 1);
    x = zeros(nC, 1);

    c1 = coordinates(elements(:, 1), :);
    d21 = coordinates(elements(:, 2), :) - c1;
    d31 = coordinates(elements(:, 3), :) - c1;

    area4 = 2 * (d21(:, 1) .* d31(:, 2) - d21(:, 2) .* d31(:, 1));

    I = reshape(elements(:, [1 2 3 1 2 3 1 2 3])', 9 * nE, 1);
    J = reshape(elements(:, [1 1 1 2 2 2 3 3 3])', 9 * nE, 1);
    a = (sum(d21 .* d31, 2) ./ area4)';
    b = (sum(d31 .* d31, 2) ./ area4)';
    c = (sum(d21 .* d21, 2) ./ area4)';
    A = [-2 * a + b + c; a - b; a - c; a - b; b; -a; a - c; -a; c];

    A = sparse(I, J, A(:), nC, nC);

    dirichlet = unique(dirichlet);

    % Check if the indices in dirichlet are within the valid range
    validDirichlet = dirichlet(dirichlet >= 1 & dirichlet <= nC);
    x(validDirichlet) = feval(uD, coordinates(validDirichlet, :));

    fsT = feval(f, c1 + (d21 + d31) / 3);
    disp('Size of elements:');
    disp(size(elements));
    disp('Size of area4:');
    disp(size(area4));
    disp('Size of fsT:');
    disp(size(fsT));
    disp('Size of x:');
    disp(size(x));
    disp('Size of A:');
    disp(size(A));

    b = accumarray(elements(:), repmat(1 / 2 * area4 .* fsT, 3, 1), [nC, 1]);
    b = [b; 0];

    if ~isempty(neumann)
        cn1 = coordinates(neumann(:, 1), :);
        cn2 = coordinates(neumann(:, 2), :);
        gmE = feval(g, (cn1 + cn2) / 2);
        b = b + accumarray(neumann(:), repmat(2 / sqrt(sum((cn2 - cn1).^2, 2)) .* gmE, 2, 1), [nC + 1, 1]);
    end

    freenodes = setdiff(1:nC, dirichlet);
    disp('Size of freenodes:');
    disp(size(freenodes));
    disp('Size of dirichlet:');
    disp(size(dirichlet));

    x(freenodes(freenodes <= nC)) = A(freenodes(freenodes <= nC), freenodes(freenodes <= nC)) \ b(freenodes(freenodes <= nC));
    x = [x; 0];

    energy = x(1:nC)' * A(1:nC, 1:nC) * x(1:nC);

end

function etaR = computeEtaR(x, coordinates, elements, dirichlet, neumann, f, g)
    [edge2nodes, element2edges, dirichlet2edges, neumann2edges] = provideGeometricData(elements, dirichlet, neumann);

    c1 = coordinates(elements(:, 1), :);
    d21 = coordinates(elements(:, 2), :) - c1;
    d31 = coordinates(elements(:, 3), :) - c1;

    area2 = d21(:, 1) .* d31(:, 2) - d21(:, 2) .* d31(:, 1);

    u21 = repmat(x(elements(:, 2)) - x(elements(:, 1)), 1, 2);
    u31 = repmat(x(elements(:, 3)) - x(elements(:, 1)), 1, 2);

    curl = (d31 .* u21 - d21 .* u31) ./ repmat(area2, 1, 2);

    dudn21 = sum(d21 .* curl, 2);
    dudn13 = -sum(d31 .* curl, 2);
    dudn32 = -(dudn13 + dudn21);

    % Initialize etaR with zeros
    etaR = zeros(size(edge2nodes, 1), 1);

    % Assign values to etaR directly, skipping zero indices
    nonZeroIndices = find(element2edges > 0);
    etaR(element2edges(nonZeroIndices)) = [dudn21; dudn32; dudn13];

    % Continue with the rest of the code
    if ~isempty(neumann)
        cn1 = coordinates(neumann(:, 1), :);
        cn2 = coordinates(neumann(:, 2), :);
        gmE = feval(g, (cn1 + cn2) / 2);
        etaR(neumann2edges) = etaR(neumann2edges) - sqrt(sum((cn2 - cn1).^2, 2)) .* gmE;
    end

    % Ensure dirichlet is a column vector
    dirichlet = dirichlet(:);

    if ~isempty(dirichlet)
        % Ensure dirichlet values are within valid range
        validDirichlet = dirichlet(dirichlet >= 1 & dirichlet <= size(coordinates, 1));

        % Check for invalid indices
        invalidIndices = setdiff(dirichlet, validDirichlet);
        if ~isempty(invalidIndices)
            disp('Invalid indices detected in dirichlet.');
            disp(['Invalid indices: ', num2str(invalidIndices')]);
            error('Please correct the indices in dirichlet.');
        end

        % Use valid indices for further computations
        dirichlet2edges(validDirichlet, validDirichlet) = 1;
    end

    % Print debugging information
    disp('etaR:');
    disp(etaR);
end

function [coordinates, elements, varargout] = refineNVB(coordinates, elements, varargin)
    if nargin < 5
        error('Insufficient input arguments.');
    end

    markedElements = varargin{end-1};
    etaR = varargin{end}; % Added to accept etaR as an input argument
    nE = size(elements, 1);

    [edge2nodes, element2edges, boundary2edges] = provideGeometricData(elements, varargin{1:end-3}, varargin{end}); % Include neumann as an additional argument

    edge2newNode = zeros(max(max(element2edges)), 1);

    if ~isempty(markedElements)
        markedEdgesIdx = element2edges(markedElements, :);
        edge2newNode(markedEdgesIdx(markedEdgesIdx >= 1)) = 1;
    end

    if any(edge2newNode)
        iterationCount = 0;

        while any(edge2newNode)
            markedEdge = zeros(size(element2edges));
            markedEdge(element2edges >= 1) = edge2newNode(element2edges(element2edges >= 1));

            if all(markedEdge(:) == 0)
                break;
            end

            swap = find(~markedEdge(:, 1) & (markedEdge(:, 2) | markedEdge(:, 3)), 1, 'first');

            if isempty(swap)
                break;
            end

             % Verify if swap corresponds to a valid element index
            if swap < 1 || swap > nE
                error('Invalid swap value computed.');
            end

            edgeIdx = element2edges(swap, 1);

            % Debugging: Output swap and edgeIdx values
            disp(['Swap value computed: ', num2str(swap)]);
            disp(['Edge index computed: ', num2str(edgeIdx)]);

            iterationCount = iterationCount + 1;

            % Check if the edge index is valid
            if edgeIdx >= 1 && edgeIdx <= numel(edge2newNode)
                if edgeIdx > numel(edge2newNode)
                    edge2newNode(edgeIdx) = numel(edge2newNode) + 1;
                else
                    edge2newNode(edgeIdx) = 1;
                end
            else
             % Handle the case when the edge index is invalid
                disp('Invalid edge index during update. Searching for a valid edge.');
                validEdge = find(edge2newNode == 0, 1);

                if isempty(validEdge)
                    disp('No valid edge found. Exiting the loop.');
                    break;
                end

                edge2newNode(validEdge) = numel(edge2newNode) + 1;
                continue;
            end

            if ~isempty(varargin{end-2}) % Check if etaR is provided
                dirichletEdges = boundary2edges{2};
                dirichletNodes = unique(edge2nodes(dirichletEdges, :));

                cn1 = coordinates(dirichletNodes, :);
                cn2 = coordinates(edge2nodes(dirichletEdges, :), :);
                gmE = feval(g, (cn1 + cn2) / 2);
                etaR(dirichletEdges) = etaR(dirichletEdges) - sqrt(sum((cn2 - cn1).^2, 2)) .* gmE;
            end
        end

        if any(edge2newNode)
            edge2newNode(edge2newNode ~= 0) = size(coordinates, 1) + (1:nnz(edge2newNode));
            idx = find(edge2newNode ~= 0);
            coordinates(edge2newNode(idx), :) = (coordinates(edge2nodes(idx, 1), :) + coordinates(edge2nodes(idx, 2), :)) / 2;

            for j = 1:nargout-2
                boundary = varargin{j};
                if ~isempty(boundary)
                    boundary2edgesIndices = [];
                    validIndices = [];
                    if ~isempty(boundary2edges{j})
                        for k = 1:numel(boundary2edges{j})
                            edgeIdx = boundary2edges{j}(k);
                            if edgeIdx >= 1 && edgeIdx <= numel(edge2newNode)
                                boundary2edgesIndices = [boundary2edgesIndices; edgeIdx];
                                validIndices = [validIndices; true];
                            else
                                error(['Invalid edge index during update: ', num2str(edgeIdx)]);
                            end
                        end
                    end

                    if any(~validIndices)
                        error('Invalid indices detected in boundary2edges{j}.');
                    end

                    newNodes = edge2newNode(boundary2edgesIndices(validIndices));

                    markedEdges = (newNodes ~= 0);
                                        if any(markedEdges)
                        validNodesIdx = find(~newNodes);
                        edge2nodes(boundary2edgesIndices(validIndices), :) = [validNodesIdx, newNodes(markedEdges)];
                    end
                end
            end

            if any(edge2newNode)
                none = ~markedEdge(:, 1);
                bisec1 = (markedEdge(:, 1) & ~markedEdge(:, 2) & ~markedEdge(:, 3));
                bisec12 = (markedEdge(:, 1) & markedEdge(:, 2) & ~markedEdge(:, 3));
                bisec13 = (markedEdge(:, 1) & ~markedEdge(:, 2) & markedEdge(:, 3));
                bisec123 = (markedEdge(:, 1) & markedEdge(:, 2) & markedEdge(:, 3));

                idx = ones(nE, 1);
                idx(bisec1) = 2;
                idx(bisec12) = 3;
                idx(bisec13) = 3;
                idx(bisec123) = 4;
                idx = [1; 1 + cumsum(idx)];

                newElements = zeros(idx(end) - 1, 3);
                newNodes = zeros(size(elements, 1), 2);  % Add this line

                newElements(idx(none), :) = elements(none, :);
                % Corrected line
                newNodes = [newNodes, zeros(size(newNodes, 1), 1)]; % Add a column of zeros to newNodes
                newNodes(edge2newNode ~= 0, 1:2) = coordinates(edge2nodes(edge2newNode ~= 0, 1), 1:2) / 2;
                newElements([idx(bisec1), 1 + idx(bisec1)], :) = [elements(bisec1, 3), elements(bisec1, 1), newNodes(bisec1, 1); elements(bisec1, 2), elements(bisec1, 3), newNodes(bisec1, 1)];
                newElements([idx(bisec12), 1 + idx(bisec12), 2 + idx(bisec12)], :) = [elements(bisec12, 3), elements(bisec12, 1), newNodes(bisec12, 1); newNodes(bisec12, 1), elements(bisec12, 2), newNodes(bisec12, 2); elements(bisec12, 3), newNodes(bisec12, 1), newNodes(bisec12, 2)];
                newElements([idx(bisec13), 1 + idx(bisec13), 2 + idx(bisec13)], :) = [newNodes(bisec13, 1), elements(bisec13, 3), newNodes(bisec13, end); elements(bisec13, 1), newNodes(bisec13, 1), newNodes(bisec13, end); elements(bisec13, 2), elements(bisec13, 3), newNodes(bisec13, 1)];
                newElements([idx(bisec123), 1 + idx(bisec123), 2 + idx(bisec123), 3 + idx(bisec123)], :) = [newNodes(bisec123, 1), elements(bisec123, 3), newNodes(bisec123, end); elements(bisec123, 1), newNodes(bisec123, 1), newNodes(bisec123, end); newNodes(bisec123, 1), elements(bisec123, 2), newNodes(bisec123, 2); elements(bisec123, 3), newNodes(bisec123, 1), newNodes(bisec123, 2)];

            else
                newElements = elements;
            end

            % Output arguments
            varargout{1} = newElements;
            varargout{2} = [];
            varargout{3} = [];
        end
    end
end

function [coordinates, elements, varargout] = coarsenNVB(N0, coordinates, elements, varargin)
    nC = size(coordinates,1);
    nE = size(elements,1);

    [I, J, nodes2edge, element2neighbours] = provideGeometricData(elements, varargin{:});

    mask = nodes2edge > 0;
    [foo{1:2}, idxIJ] = find(nodes2edge);
    [foo{1:2}, neighbourIJ] = find(mask + mask.*sparse(J,I,[1:nE,1:nE,1:nE]'));
    element2neighbours(idxIJ) = neighbourIJ - 1;
    element2neighbours = reshape(element2neighbours,nE,3);

    marked = zeros(nE,1);
    marked(varargin{end}) = 1;
    newestNode = unique(elements((marked & elements(:,3) > N0),3));
    valence = accumarray(elements(:),1,[nC 1]);
    markedNodes = zeros(nC,1);
    markedNodes(newestNode((valence(newestNode) == 2 | valence(newestNode) == 4))) = 1;

    idx = find(markedNodes(elements(:,3)) & (element2neighbours(:,3) > (1:nE)'))';
    markedElements = zeros(nE,1);
    markedElements(idx) = 1;

    for element = idx
        if markedElements(element)
            markedElements(element2neighbours(element,3)) = 0;
        end
    end

    idx = find(markedElements);
    brother = element2neighbours(idx,3);
    elements(idx,[1 3 2]) = [elements(idx,[2 1]) elements(brother,1)];

    activeNodes = find(~markedNodes);
    coordinates = coordinates(activeNodes,:);

    coordinates2newCoordinates = zeros(1,nC);
    coordinates2newCoordinates(activeNodes) = 1:length(activeNodes);

    elements(brother,:) = [];
    elements = coordinates2newCoordinates(elements);

    for j = 1:nargout-2
        boundary = varargin{j};
        if ~isempty(boundary)
            node2boundary = zeros(nC,2);
            node2boundary(boundary(:,1),1) = 1:size(boundary,1);
            node2boundary(boundary(:,2),2) = 1:size(boundary,1);
            idx = (markedNodes & node2boundary(:,2));
            boundary(node2boundary(idx,2),2) = boundary(node2boundary(idx,1),2);
            boundary(node2boundary(idx,1),2) = 0;
            varargout{j} = coordinates2newCoordinates(boundary(find(boundary(:,2)),:));
        else
            varargout{j} = [];
        end
    end
end

function [x,coordinates,elements,indicators] = adaptiveAlgorithm(coordinates,elements,dirichlet,neumann,f,g,uD,nEmax,theta)
    while 1
        x = solveLaplace(coordinates,elements,dirichlet,neumann,f,g,uD);
        indicators = computeEtaR(x,coordinates,elements,dirichlet,neumann,f,g);

        if size(elements,1) >= nEmax
            break
        end

        [indicators,idx] = sort(indicators,'descend');
        sumeta = cumsum(indicators);
        ell = find(sumeta >= sumeta(end) * theta,1);
        marked = idx(1:ell);

        [coordinates,elements,dirichlet,neumann] = refineNVB(coordinates,elements,dirichlet,neumann,marked);
    end
end

function etaC = findEtaCoarse(elements, marked)
    etaC = zeros(size(elements, 1), 1);
    etaC(marked) = 1;
end

function delta_t = computeDeltaT(x, coordinates, elements, dirichlet, neumann, f, g, theta, nEmax)
    indicators = computeEtaR(x, coordinates, elements, dirichlet, neumann, f, g);

    if size(elements, 1) >= nEmax
        delta_t = inf;
        return;
    end

    [indicators, idx] = sort(indicators, 'descend');
    sumeta = cumsum(indicators);
    ell = find(sumeta >= sumeta(end) * theta, 1);

    if isempty(ell)
        delta_t = inf;
    else
        marked = idx(1:ell);
        delta_t = 1 / sqrt(sum(1./indicators(marked)));
    end
end
