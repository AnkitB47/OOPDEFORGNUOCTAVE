classdef SlicedSquareModified < grid2D
    properties

        outer_boundary % Outer boundary points
        inflow_boundary % Inflow boundary points
        outflow_boundary % Outflow boundary points

    end

    methods(Access = public)
        function obj = SlicedSquareModified()
            fd = @(p) ddiff(ddiff(drectangle(p,-1,1,-1,1),...
                drectangle(p,-0.55,-0.45,-1,0)),...
                drectangle(p,0.45,0.55,-1,0));
            fixPoints1 = [linspace(-1,-0.55,15);-ones(1,15)]';
            fixPoints2 = [linspace(-0.45,0.45,20);-ones(1,20)]';
            fixPoints3 = [linspace(0.55,1,15);-ones(1,15)]';

            fh = @(p) 0.3+drectangle(p,-0.55,-0.45,-1,0)+drectangle(p,0.45,0.55,-1,0);
            [p,e,t] = distmesh2d(fd,fh,.05,...
                                [-1,-1;1,1],...
                                [-1,1;...
                                  1,1;...
                                  fixPoints1;...
                                  fixPoints2;...
                                  fixPoints3;...
                                -0.55,0;...
                                -0.45,0;...
                                0.45,0;...
                                0.55,0]);

            obj.p = p';
            obj.t = [t';ones(1,size(t',2))];
            obj.e = e';

            % Identify outer boundary using convex hull
            obj.identifyOuterBoundary();
        end

        function identifyBoundaries(obj)
            % Convert mesh points to column vectors
            points = obj.p';

            % Find outer boundary using convex hull modified
            obj.outer_boundary = convexHullModified(points);

            % Check if outer boundary is empty
            if isempty(obj.outer_boundary)
                error('Outer boundary is empty. Check the geometry definition or the convexHullModified function.');
            end

            % Initialize arrays for inflow, outflow, and boundary points
            obj.inflow_boundary = [];
            obj.outflow_boundary = [];

            % Define tolerance for considering points as inflow or outflow
            inflow_tolerance = 1e-3;
            outflow_tolerance = 1 - inflow_tolerance;

            % Iterate through all mesh points
            for i = 1:size(points, 1)
                point = points(i, :);

                % Check if point is near any of the four segments of the outer boundary
                dist1 = pointToLineSegmentDistance(point, obj.outer_boundary(1,:), obj.outer_boundary(2,:));
                dist2 = pointToLineSegmentDistance(point, obj.outer_boundary(2,:), obj.outer_boundary(3,:));
                dist3 = pointToLineSegmentDistance(point, obj.outer_boundary(3,:), obj.outer_boundary(4,:));
                dist4 = pointToLineSegmentDistance(point, obj.outer_boundary(4,:), obj.outer_boundary(1,:));

                % Determine classification based on distance to segments
                if dist1 < inflow_tolerance || dist2 < inflow_tolerance || dist3 < inflow_tolerance || dist4 < inflow_tolerance
                    obj.inflow_boundary = [obj.inflow_boundary; point];
                elseif dist1 > outflow_tolerance && dist2 > outflow_tolerance && dist3 > outflow_tolerance && dist4 > outflow_tolerance
                    obj.outflow_boundary = [obj.outflow_boundary; point];
                end
            end
        end

        function dist = pointToLineSegmentDistance(point, line_start, line_end)
            % Vector from line start to point
            v1 = point - line_start;

            % Vector from line start to line end
            v2 = line_end - line_start;

            % Length of v2
            len_v2 = norm(v2);

            % Project v1 onto v2
            proj = dot(v1, v2) / len_v2;

            % If the projection is less than 0, the closest point is the line start
            if proj < 0
                dist = norm(point - line_start);
            % If the projection is greater than the length of v2, the closest point is the line end
            elseif proj > len_v2
                dist = norm(point - line_end);
            % Otherwise, calculate the distance to the line
            else
                % Closest point on the line
                closest_point = line_start + (proj / len_v2) * v2;
                % Distance from the point to the closest point on the line
                dist = norm(point - closest_point);
            end
        end

        function plotHighlightedMesh(obj)
            % Plot the mesh
            obj.plot();

            % Highlight inflow, outflow, and boundary regions
            hold on;
            scatter(obj.p(obj.inflow_indices, 1), obj.p(obj.inflow_indices, 2), 'g', 'filled');
            scatter(obj.p(obj.outflow_indices, 1), obj.p(obj.outflow_indices, 2), 'r', 'filled');
            scatter(obj.p(obj.boundary_indices, 1), obj.p(obj.boundary_indices, 2), 'b', 'filled');
            hold off;

            % Add legend
            legend('Mesh', 'Inflow', 'Outflow', 'Boundary');
        end

    end
end

