% Define the distance function for a quadrant circle
fd_circle = @(p, varargin) sqrt(p(:, 1).^2 + p(:, 2).^2) - 1;

% Define the edge length function for uniform mesh
fh_non_uniform =  @(p) 0.05 + 5 * abs((p(:,1).^2 - p(:,2).^2));  % Set desired edge length

% Example values, adjust as needed
geps = 1e-5;
h0 = 0.1;
bbox = [-1, -1; 1, 1];
pfix = [];  % Your fixed points array, if any
max_iterations = 500;

% Set bounding box
bbox = [-1, -1; 1, 1];

% No fixed node positions
pfix = [];

% Additional parameters (empty for now)
varargin = [];

% Create initial mesh using distmesh
[p_initial, t_initial, coordinates, elements] = do_distmesh2d_polyhedral(fd_circle, fh_non_uniform, h0, bbox, pfix, geps, max_iterations, varargin);

% Assign the mesh nodes obtained from do_distmesh2d_polyhedral to N0
N0 = size(p_initial,1);

% Define the source term function handle
f = @(p) ones(size(p, 1), 1); % For example, a constant source term

% Define the Dirichlet boundary condition function handle
uD = @(p) zeros(size(p, 1), 1); % For example, zero Dirichlet boundary condition

% Define the Neumann boundary condition function handle
g = @(p) zeros(size(p, 1), 1); % For example, zero Neumann boundary condition

theta = 0.5; % Dörfler criterion parameter

% Assuming 'elements' contains the mesh elements
boundaries = computeBoundaries(elements);

% Call adaptiveMeshRef with all required function handles
[refined_coordinates, refined_elements] = adaptiveMeshRef(N0, coordinates, elements, [], [], uD, boundaries, f, g);

% Plot the initial and refined meshes side by side
figure;

% Plot the initial mesh
subplot(1, 2, 1);
trimesh(t_initial, p_initial(:,1), p_initial(:,2));
title('Initial Mesh');

% Plot the refined mesh
subplot(1, 2, 2);
trimesh(refined_elements, refined_coordinates(:,1), refined_coordinates(:,2));
title('Refined Mesh');

