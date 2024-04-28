classdef VibrationIsolator < grid2D
    methods(Access = public)
        function obj =  VibrationIsolator(h)
            if nargin == 0
                h = 0.1;
            end
            h = min(max(h,0.025),0.2);
            s = 0:h:2*pi ;
            s2 = 0:8*h:2*pi;   
            s3 = 0:3*h:2*pi;
            obj.freeGeometry([sin(s);-cos(s)],...
                [0.125*sin(s2)+0.75;0.125*cos(s2)],...
                [0.125*sin(s2)-0.75;0.125*cos(s2)],...
                [0.125*sin(s2);0.125*cos(s2)+0.75],...
                [0.125*sin(s2);0.125*cos(s2)-0.75],...
                [0.2*sin(s3);0.2*cos(s3)]);
            obj.jiggleMesh;
        end
    end
end
