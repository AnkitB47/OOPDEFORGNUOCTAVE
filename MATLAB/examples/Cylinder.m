classdef Cylinder < grid3D    
% Cylinder constructor for cylinder class
%
% Creates a cylinder geometrie using tetrahedrons.
% Arguments are R,H,[h][,x,y,z]
%
    
    %%
% 
    
    methods(Access = public)
        function obj = Cylinder(varargin)
            obj = obj@grid3D(); 
            obj.cylinder(varargin{:});
        end
    end
    
end

