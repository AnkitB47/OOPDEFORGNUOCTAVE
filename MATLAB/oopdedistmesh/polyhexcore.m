classdef polyhexcore
    % Copyright (C) 2024 Ankit Bhardwaj
    properties
        BoundingBox
        FixedPoints
    end

    methods
        function obj = polyhexcore(dim)
            % Constructor
            if nargin < 1
                dim = 3;
            end

            % Set up the bounding box
            if dim == 2
                obj.BoundingBox = [-1, -1; 1, 1];
            elseif dim == 3
                obj.BoundingBox = [-1, -1, -1; 1, 1, 1];
            else
                error('Invalid dimension. Use 2 or 3.');
            end

            % Set up the fixed points
            if dim == 2
                obj.FixedPoints = [];
            elseif dim == 3
                obj.FixedPoints = [0, 0, 0];
            end
        end

        function generateMesh(obj)
            % Generate mesh using distmesh2d
            %fd = @(p) dsphere(p, 0, 0, 0.5);

            % Call distmesh2d with the appropriate distance function
            if size(obj.BoundingBox, 2) == 2
                fd = @(p) dcircle(p, 0, 0, 0.5);
                [p, ~, t] = distmesh2d(fd, @huniform, 0.1, obj.BoundingBox, obj.FixedPoints);
                trimesh(t, p(:, 1), p(:, 2));
                axis equal;
            elseif size(obj.BoundingBox, 2) == 3
                fd = @(p) dsphere(p, 0, 0, 0, 0.5);
                [p, t] = distmeshnd(fd, @huniform, 0.1, obj.BoundingBox, obj.FixedPoints);
                tetramesh(t, p);
                axis equal;
            end

            % Plot the mesh
            %if size(obj.BoundingBox, 2) == 2
            %    trimesh(t, p(:, 1), p(:, 2));
            %elseif size(obj.BoundingBox, 2) == 3
            %    tetramesh(t, p);
            %end
            %axis equal;
        end
    end
end

