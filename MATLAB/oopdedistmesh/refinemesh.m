function [x, coordinates, elements, indicators] = adaptiveMeshRefinement(coordinates, elements, dirichlet, neumann, f, g, uD, tmax, nEmax, theta, tau, sigma)
    n = 0;
    t = 0;
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

            % Compute discrete solution on coarser mesh
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

% Define other functions used in the algorithm (solveLaplace, computeEtaR, refineNVB, coarsenNVB, findEtaCoarse) based on your provided code snippets.
function [x,energy] = solveLaplace(coordinates,elements,dirichlet,neumann,f,g,uD)
    nE = size(elements,1);
    nC = size(coordinates,1);
    x = zeros(nC,1);

    c1 = coordinates(elements(:,1),:);
    d21 = coordinates(elements(:,2),:) - c1;
    d31 = coordinates(elements(:,3),:) - c1;

    area4 = 2 * (d21(:,1) .* d31(:,2) - d21(:,2) .* d31(:,1));

    I = reshape(elements(:, [1 2 3 1 2 3 1 2 3])', 9 * nE, 1);
    J = reshape(elements(:, [1 1 1 2 2 2 3 3 3])', 9 * nE, 1);
    a = (sum(d21 .* d31,2) ./ area4)';
    b = (sum(d31 .* d31,2) ./ area4)';
    c = (sum(d21 .* d21,2) ./ area4)';
    A = [-2*a+b+c;a-b;a-c;a-b;b;-a;a-c;-a;c];

    A = sparse(I, J, A(:));

    dirichlet = unique(dirichlet);
    x(dirichlet) = feval(uD, coordinates(dirichlet,:));

    fsT = feval(f, c1 + (d21 + d31) / 3);
    b = accumarray(elements(:), repmat(1 / 2 * area4 .* fsT, 3, 1), [nC 1]) - A * x;

    if ~isempty(neumann)
        cn1 = coordinates(neumann(:,1),:);
        cn2 = coordinates(neumann(:,2),:);
        gmE = feval(g, (cn1 + cn2) / 2);
        b = b + accumarray(neumann(:), repmat(2 / sqrt(sum((cn2 - cn1).^2,2)) .* gmE, 2, 1), [nC 1]);
    end

    freenodes = setdiff(1:nC, dirichlet);
    x(freenodes) = A(freenodes, freenodes) \ b(freenodes);

    energy = x' * A * x;
end

function etaR = computeEtaR(x,coordinates,elements,dirichlet,neumann,f,g)
    [edge2nodes,element2edges,dirichlet2edges,neumann2edges] = provideGeometricData(elements, dirichlet, neumann);

    c1 = coordinates(elements(:,1),:);
    d21 = coordinates(elements(:,2),:) - c1;
    d31 = coordinates(elements(:,3),:) - c1;

    area2 = d21(:,1) .* d31(:,2) - d21(:,2) .* d31(:,1);

    u21 = repmat(x(elements(:,2)) - x(elements(:,1)), 1, 2);
    u31 = repmat(x(elements(:,3)) - x(elements(:,1)), 1, 2);
    curl = (d31 .* u21 - d21 .* u31) ./ repmat(area2, 1, 2);

    dudn21 = sum(d21 .* curl, 2);
    dudn13 = -sum(d31 .* curl, 2);
    dudn32 = -(dudn13 + dudn21);
    etaR = accumarray(element2edges(:), [dudn21;dudn32;dudn13], [size(edge2nodes,1) 1]);

    if ~isempty(neumann)
        cn1 = coordinates(neumann(:,1),:);
        cn2 = coordinates(neumann(:,2),:);
        gmE = feval(g, (cn1 + cn2) / 2);
        etaR(neumann2edges) = etaR(neumann2edges) - sqrt(sum((cn2 - cn1).^2,2)) .* gmE;
    end

    etaR(dirichlet2edges) = 0;
    etaR = sum(etaR(element2edges).^2, 2) + 0.5 * area2 .* feval(f, c1 + (d21 + d31) / 3).^2;
end

function [coordinates, newElements, varargout] = refineNVB(coordinates, elements, varargin)
    markedElements = varargin{end};
    nE = size(elements,1);

    [edge2nodes, element2edges, boundary2edges] = provideGeometricData(elements, varargin{1:end-1});

    edge2newNode = zeros(max(max(element2edges)), 1);
    edge2newNode(element2edges(markedElements,:)) = 1;

    swap = 1;
    while ~isempty(swap)
        markedEdge = edge2newNode(element2edges);
        swap = find(~markedEdge(:,1) & (markedEdge(:,2) | markedEdge(:,3)));
        edge2newNode(element2edges(swap,1)) = 1;
    end

    edge2newNode(edge2newNode ~= 0) = size(coordinates,1) + (1:nnz(edge2newNode));
    idx = find(edge2newNode);
    coordinates(edge2newNode(idx),:) = (coordinates(edge2nodes(idx,1),:) + coordinates(edge2nodes(idx,2),:)) / 2;

    for j = 1:nargout-2
        boundary = varargin{j};
        if ~isempty(boundary)
            newNodes = edge2newNode(boundary2edges{j});
            markedEdges = find(newNodes);
            if ~isempty(markedEdges)
                boundary = [boundary(~newNodes,:); boundary(markedEdges,1), newNodes(markedEdges); newNodes(markedEdges), boundary(markedEdges,2)];
            end
        end
        varargout{j} = boundary;
    end

    newNodes = edge2newNode(element2edges);
    markedEdges = (newNodes ~= 0);
    none = ~markedEdges(:,1);
    bisec1 = (markedEdges(:,1) & ~markedEdges(:,2) & ~markedEdges(:,3));
    bisec12 = (markedEdges(:,1) & markedEdges(:,2) & ~markedEdges(:,3));
    bisec13 = (markedEdges(:,1) & ~markedEdges(:,2) & markedEdges(:,3));
    bisec123 = (markedEdges(:,1) & markedEdges(:,2) & markedEdges(:,3));

    idx = ones(nE,1);
    idx(bisec1) = 2;
    idx(bisec12) = 3;
    idx(bisec13) = 3;
    idx(bisec123) = 4;
    idx = [1; 1 + cumsum(idx)];

    newElements = zeros(idx(end) - 1, 3);
    newElements(idx(none),:) = elements(none,:);
    newElements([idx(bisec1),1 + idx(bisec1)],:) = [elements(bisec1,3), elements(bisec1,1), newNodes(bisec1,1); elements(bisec1,2), elements(bisec1,3), newNodes(bisec1,1)];
    newElements([idx(bisec12),1 + idx(bisec12),2 + idx(bisec12)],:) = [elements(bisec12,3), elements(bisec12,1), newNodes(bisec12,1); newNodes(bisec12,1), elements(bisec12,2), newNodes(bisec12,2); elements(bisec12,3), newNodes(bisec12,1), newNodes(bisec12,2)];
    newElements([idx(bisec13),1 + idx(bisec13),2 + idx(bisec13)],:) = [newNodes(bisec13,1), elements(bisec13,3), newNodes(bisec13,3); elements(bisec13,1), newNodes(bisec13,1), newNodes(bisec13,3); elements(bisec13,2), elements(bisec13,3), newNodes(bisec13,1)];
    newElements([idx(bisec123),1 + idx(bisec123),2 + idx(bisec123),3 + idx(bisec123)],:) = [newNodes(bisec123,1), elements(bisec123,3), newNodes(bisec123,3); elements(bisec123,1), newNodes(bisec123,1), newNodes(bisec123,3); newNodes(bisec123,1), elements(bisec123,2), newNodes(bisec123,2); elements(bisec123,3), newNodes(bisec123,1), newNodes(bisec123,2)];
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
