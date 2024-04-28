classdef channelWithCavity < grid2D
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

     

    methods
        function obj = channelWithCavity(h)
            if nargin == 0
                h = 0.1;
            end
            x = [0 0 1 2 2 3 3 4 4 3 2 1];
            y = [2 1 1 1 0 0 1 1 2 2 2 2];
            obj.freeGeometry([x;y]);
             
            obj.e(end,:) = [1 2 2 2 2 2 2 3 2 2 2 2]
            while h < obj.hmax
                obj.refineMesh;
            end
        end
             
             
    end
         
     
end