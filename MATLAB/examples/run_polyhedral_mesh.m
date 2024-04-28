% Specify the dimensions of the rectangle
width = 2;
height = 1;

% Number of nodes along each dimension
M = 10;

% Generate vertices for the rectangle
vertices = [-width / 2, -height / 2;
             width / 2, -height / 2;
             width / 2, height / 2;
            -width / 2, height / 2];

% Set N proportional to the number of nodes
N = 4 * M;

% Create edges and faces for the rectangle
[edges, faces] = computeEdgesAndFaces(vertices);

% Set the desired element size for the mesh
h0 = 0.1;

% Define the fd function for the rectangle
fd = @(p) fd_rectangle(p, -width / 2, width / 2, -height / 2, height / 2);

% Define the fh function (uniform mesh size)
fh = @(p) h0 * ones(size(p, 1), 1);

% Call distmesh2d_polyhedral with edges and faces, and pass fd and fh
[p, e, t] = distmesh2d_polyhedral(fd, fh, h0, vertices, edges, faces);

% Check if t is defined before attempting to plot the mesh
if exist('t', 'var') && ~isempty(t)
    % Plot the mesh
    figure;
    triplot(t, p(:, 1), p(:, 2), 'k-');
    axis equal;
    xlabel('X');
    ylabel('Y');
    title('Polyhedral Mesh for Rectangle');
else
    disp('Mesh generation failed due to too few points.');
end

