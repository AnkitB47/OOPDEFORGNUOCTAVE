classdef PolygonGeometry < grid2D
    methods
        function obj = PolygonGeometry(polygon,h)
            obj.polygonGeometry(polygon,h);


            % Correct boundary segment numbering
            [m,n] = size(polygon)
            % force polygon to be a row matrix.
            if m > n
                polygon = polygon';
            end
            % Identify corner points from polygon
            cornerPoints =  zeros(1,size(polygon,2));
            for k = 1:max(size(polygon))
                n1 = obj.p(1,:)==polygon(1,k);
                n2 = obj.p(2,:)==polygon(2,k);
                cornerPoints(k) = find(n1.*n2 == 1);
            end
            % Keep number of corner points and extend cornerpoint by the
            % first one to wrap arround.
            NCP = length(cornerPoints);
            cornerPoints = [cornerPoints cornerPoints(1)];
            % We start with boundary segment 1
            BS = 1;
            % For number of segments in polygon we look for the elements
            % belonging to one segment
            for k1 = 1:NCP
                akt = cornerPoints(k1)
                while cornerPoints(k1+1)~=akt
                    k2 = find(akt==obj.e(1,:));
                    obj.e(end,k2) = BS;
                    akt = obj.e(2,k2);
                end
                BS = BS + 1;
            end
        end
    end
end

