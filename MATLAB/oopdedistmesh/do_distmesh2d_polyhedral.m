function [p, t, coordinates, elements] = do_distmesh2d_polyhedral(fd, fh, h0, bbox, pfix, geps, max_iterations, varargin)
    Fscale = 1.0;
    deltat = 0.2;
    dptol = 0.01;
    ttol = 0.1;
    % Step 1: Create uniform distribution of nodes
    [x, y] = meshgrid(bbox(1,1):h0:bbox(2,1), bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
    x(2:2:end,:) = x(2:2:end,:) + h0/2; % Shift even rows
    p = [x(:), y(:)]; % List of node coordinates

    % Step 2: Remove nodes outside the desired geometry
    d = feval(fd, p);
    p = p(all(d < -geps, 2), :);

    % Rejection method
    if nargin > 4 && ~isempty(pfix)
        p = [pfix; p];
    end

    % Display sizes after modifications
    disp('Size of p:');
    disp(size(p));
    disp('Size of d:');
    disp(size(d));

    N = size(p, 1);

    % Initialize pold for the first iteration
    pold = inf;

    for count = 1:max_iterations
        fprintf('Iteration: %d of %d\n', count, max_iterations);
        fprintf('Number of points: %d\n', size(p, 1));

        % Step 3: Retriangulation by the Delaunay algorithm
        if count == 1 || max(sqrt(sum((p - pold).^2, 2))/h0) > ttol
            pold = p; % Save current positions

            % Triangulation
            try
                t = delaunayn(p);
                pmid = (p(t(:,1),:) + p(t(:,2),:) + p(t(:,3),:))/3;
                t = t(feval(fd, pmid) < -geps, :);
            catch
                fprintf('Error in Delaunay triangulation\n');
                break;
            end

            if isempty(t)
                fprintf('Error: No valid triangles in the mesh\n');
                break;
            end

            % Assign t to elements
            elements = t;

            % Step 4: Create list of edges (bars)
            bars = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])]; % Interior bars duplicated
            bars = unique(sort(bars, 2), 'rows'); % Bars as node pairs

            % Step 5: Graphical output (optional)
            % trimesh(t, p(:,1), p(:,2), zeros(N, 1));
            % view(2), axis equal, axis off, drawnow;

            % Step 6: Compute bar lengths and forces
            barvec = p(bars(:,1),:) - p(bars(:,2),:);
            L = sqrt(sum(barvec.^2, 2));
            hbars = feval(fh, (p(bars(:,1),:) + p(bars(:,2),:))/2);
            L0 = hbars * Fscale * sqrt(sum(L.^2) / sum(hbars.^2));
            F = max(L0 - L, 0);

            % Step 7: Move mesh points based on bar forces
            Fvec = F./L .* [ones(size(F)), ones(size(F))] .* barvec; % Bar forces (x,y components)
            Ftot = full(sparse(bars(:,[1,1,2,2]), ones(size(F))*[1,2,1,2], [Fvec, -Fvec], N, 2));
            Ftot(1:size(pfix, 1), :) = 0; % Force = 0 at fixed points
            p = p + deltat * Ftot;

            % Step 8: Project points back to the boundary if outside
            d = feval(fd, p);
            ix = d > 0; % Find points outside (d > 0)
            deps = 1e-6;
            dgradx = (feval(fd, [p(ix,1)+deps, p(ix,2)], varargin{:}) - d(ix))/deps; % Numerical gradient
            dgrady = (feval(fd, [p(ix,1), p(ix,2)+deps], varargin{:}) - d(ix))/deps; % Numerical gradient
            p(ix,:) = p(ix,:) - [d(ix).*dgradx, d(ix).*dgrady]; % Project back to boundary

            % Step 9: Termination criterion
            if max(sqrt(sum(deltat * Ftot(d<-geps,:).^2, 2))/h0) < dptol
                fprintf('Mesh has converged. Exiting loop.\n');
                break;
            end
        end
    end
    coordinates = p;
    elements = t;
    fprintf('Termination conditions met. Exiting loop.\n');
    fprintf('Mesh generation completed.\n');
end
