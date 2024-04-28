% Define grid parameters
nElements = 10; % Number of elements
nPoints = nElements + 1; % Number of points
x = linspace(0, 1, nPoints); % Grid points
c = ones(nElements, 1); % Convection coefficients

% Create instance of Lagrange11D class
lagrange = Lagrange11D();

% Compute gradient matrices
DX = lagrange.gradientMatrices(x);

% Compute flux through element edges
flux_edges = Lagrange11D.fluxThroughEdges(lagrange, x, c);

% Compute flux jumps
jumps = Lagrange11D.fluxJumps(lagrange, flux_edges, []);

% Display results or perform further computations
disp('Gradient Matrices:');
disp(DX);
disp('Flux Through Element Edges:');
disp(flux_edges);
disp('Flux Jumps:');
disp(jumps);

