classdef Triangle < grid2D
    methods(Access = public)
        function obj = Triangle(hmax)
            
            switch nargin
                case 0
                case 1
                    obj.p = [0 1 0 0.5 0.5 0
                             0 0 1 0   0.5 0.5];

                    obj.e = [ 1   4   2   5   3   6
                              4   2   5   3   6   1
                              0   0.5 0   0.5 0 0.5
                              0.5 1   0.5 1   0.5 1
                              1   1   2   2   3   3];

                    obj.t = [1 2 4 5
                             4 5 5 3
                             6 4 6 6
                             1 1 1 1];
                         
                 while max(obj.triangleDiameters)>hmax
                     obj.refineMesh;
                 end
            end
        end
    end
end