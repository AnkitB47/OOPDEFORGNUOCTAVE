function [coordinates, elements, t] = d_polyhedral(fd, fh, h0, bbox, pfix, max_iters, varargin)
    % Initialize variables
    coordinates = [];
    elements = [];

    % Other parameters
    geps = 1e-5;

    % DistMesh algorithm
    [p, t] = do_distmesh2d_polyhedral(fd, fh, h0, bbox, pfix, varargin{:});

    % Clean-up unused nodes and elements
    [coordinates, elements] = cleanup_mesh(p, t);
end
