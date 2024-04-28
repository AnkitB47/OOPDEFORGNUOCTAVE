classdef Part3D < grid3Dpr
    methods(Access = public)
        function obj = Part3D(h)
            if nargin == 0
                h = 10;
            end
            obj = obj@grid3Dpr;
            h = min(max(h,3),30);
            
            g2D = Part(h);
            h = max(g2D.triangleDiameters);
            n = ceil(10/h)+1;
            obj.extrude(g2D,linspace(0,10,n));
        end
    end
end
    