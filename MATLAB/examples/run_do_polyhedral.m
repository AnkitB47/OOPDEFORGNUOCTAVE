% Example values, adjust as needed
geps = 1e-5;
h0 = 0.1;
bbox = [0, 0; 1, 1];
pfix = [];  % Your fixed points array, if any
max_iterations = 100;


% Define the distance function for a quadrant circle
fd_circle = @(p, varargin) sqrt(p(:,1).^2 + p(:,2).^2) - 1;

% Define the edge length function for uniform mesh
fh_uniform = @(p) 0.05 * ones(size(p, 1), 1); % Set desired edge length

% Define the size function for a non-uniform triangular mesh
fh_non_uniform = @(p) 0.05 + 5 * abs((p(:,1).^2 - p(:,2).^2)); % Adjust based on distance from circle boundary

% Set initial edge length
h0 = 0.1;

% Set bounding box
bbox = [-1, -1; 1, 1];

% No fixed node positions
pfix = [];

% Additional parameters (empty for now)
varargin = [];

% Call the meshing function for non-uniform mesh
[p_uniform, t_uniform] = do_distmesh2d_polyhedral(fd_circle, fh_uniform, h0, bbox, pfix, geps, max_iterations, varargin);
[p_non_uniform, t_non_uniform] = do_distmesh2d_polyhedral(fd_circle, fh_non_uniform, h0, bbox, pfix, geps, max_iterations, varargin);

% Plot the uniform mesh
figure;
subplot(1, 2, 1);
trimesh(t_uniform, p_uniform(:,1), p_uniform(:,2));
title('Uniform Mesh');

% Plot the non-uniform mesh
subplot(1, 2, 2);
trimesh(t_non_uniform, p_non_uniform(:,1), p_non_uniform(:,2));
title('Non-uniform Mesh');

