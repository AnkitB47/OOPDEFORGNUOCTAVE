classdef Cavity < grid2D
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

     

    methods
        function obj = Cavity
            P = [0 1 2 2 1 0
                0 0 0 1 1 1];
            
             obj.freeGeometry(P);
             obj.e(end,:) = [1 1 1 1 2 1];
        end

         
    end
end