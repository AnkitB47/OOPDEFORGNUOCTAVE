classdef Star < grid2D
    methods(Access = public)
        function obj = Star(h)
            if nargin == 0
                h = 0.1;
            end
            h = min(max(h,0.0125),0.2);
            k = 5;    
            r = 0.3;      
            s = 0:h:2*pi;
            obj.freeGeometry(...
                [r*(k-1)*cos(s)+r*cos((k-1)*s);...
                r*(k-1)*sin(s)-r*sin((k-1)*s)],...
                [[linspace(-0.3,0.25,10),0.3*ones(1,10),...
                linspace(0.3,-0.25,10),-0.3*ones(1,10)];...
                [-0.3*ones(1,10),linspace(-0.3,0.25,10),...
                0.3*ones(1,10),linspace(0.3,-0.25,10)]]);
            obj.jiggleMesh;
        end
    end
end
