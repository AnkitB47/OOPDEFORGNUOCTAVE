classdef Sheet < grid3Dpr
    methods(Access = public)
        function obj = Sheet(h)
            if nargin == 0
                h = 0.00125;
            end
            obj = obj@grid3Dpr;
            
            
            g2D = Rectangle(0,2,0,0.005,h);
           
             
            obj.extrude(g2D,[0 0.02 0.09 0.195  0.28 0.3]);
        end
    end
end
    