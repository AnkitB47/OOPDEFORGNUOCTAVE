classdef Oval < grid2D
    methods(Access = public)
        function obj =  Oval(h)
            if nargin == 0
                h = 0.2;
            end
            h = min(max(h,0.1),0.3);
            s = 0:h:2*pi;            
            obj.freeGeometry([sin(s)+0.2*sin(2*s);-cos(s)],h);
            obj.jiggleMesh;
        end 
    end
end
