classdef CylinderP < grid3Dpr    
% Cylinder constructor for cylinder class
%
% Creates a cylinder geometrie using prism elements.
% Arguments are R,H,[h][,x,y,z]
%
    
%%
% 
    
    methods(Access = public)
        function obj = CylinderP(varargin)
            obj = obj@grid3Dpr(); 
            obj.cylinder(varargin{:});
        end
    end
    
end

