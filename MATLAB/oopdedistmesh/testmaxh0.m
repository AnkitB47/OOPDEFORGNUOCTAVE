% Define the distance function for a quadrant circle
fd_circle = @(p, varargin) sqrt(p(:, 1).^2 + p(:, 2).^2) - 1;

% Define the edge length function for uniform mesh
fh_non_uniform =  @(p) 0.05 + 5 * abs((p(:,1).^2 - p(:,2).^2));  % Set desired edge length

% Example values
geps = 1e-5;
h0 = 0.11;
bbox = [-1, -1; 1, 1];
pfix = [];  % Your fixed points array, if any
max_iterations = 500;

% Define the source term function handle
f = @(p) ones(size(p, 1), 1); % For example, a constant source term

% Define the Dirichlet boundary condition function handle
uD = @(p) zeros(size(p, 1), 1); % zero Dirichlet boundary condition

% Define the Neumann boundary condition function handle
g = @(p) zeros(size(p, 1), 1); % zero Neumann boundary condition

theta = 0.5; % Dörfler criterion parameter

[p_initial, ~, coordinates, elements] = do_distmesh2d_polyhedral(fd_circle, fh_non_uniform, h0, bbox, pfix, geps, max_iterations, []);

% Call the function to find max allowable h0
max_h0 = findMaxAllowableH0(h0, fd_circle, fh_non_uniform, bbox, pfix, geps, max_iterations, f, g, uD, size(p_initial,1), theta);

disp(['Maximum allowable h0: ', num2str(max_h0)]);

