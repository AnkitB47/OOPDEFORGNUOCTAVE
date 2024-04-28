classdef Part < grid2D
    %PART class that defines a construction part   
    
    methods(Access = public)        
        function obj = Part(h)
            if nargin == 0
               %
            end
            h = min(max(h,3),30);
             
            N1 = 20; N2= 8; N3 = 12; N4 = 12; N5 = 8; N6 = 20; 
            phi = linspace(0,2*pi,20);
            x1 = linspace(0,100,N1); x1(end) = [];
            x2 = 100*ones(1,N2); x2(end) = [];
            x3 = linspace(100,40,N3); x3(end) = [];
            x4 = 40*ones(1,N4); x4(end) = [];
            x5 = linspace(40,0,N5); x5(end) = [];
            x6 = zeros(1,N6); x6(end) = [];
            y1 = zeros(1,N1); y1(end) = [];
            y2 = linspace(0,40,N2); y2(end) = [];
            y3 = 40*ones(1,N3); y3(end) = [];
            y4 = linspace(40,100,N4); y4(end) = [];
            y5 = 100*ones(1,N5); y5(end) = [];
            y6 = linspace(100,0,N6); y6(end) = [];
             
            obj.freeGeometry([x1 x2 x3 x4 x5 x6;...
                y1 y2 y3 y4 y5 y6],...
                [8*sin(phi)+50;8*cos(phi)+20],...
                [8*sin(phi)+80;8*cos(phi)+20],...
                [8*sin(phi)+20;8*cos(phi)+50],...
                [8*sin(phi)+20;8*cos(phi)+80]);
            % only one boundary segment
            obj.e(5,:) =  ones(size(obj.e(5,:)));  
            
            while max(obj.triangleDiameters) > h 
                obj.refineMesh;
            end
            obj.jiggleMesh;
        end
    end         
end
