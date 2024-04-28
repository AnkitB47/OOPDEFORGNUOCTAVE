classdef Cube < grid3D    
     
    
    methods(Access = public)
        function obj = Cube(l,hmax)
            obj = obj@grid3D(); 
            obj.bar(0,l,0,l,0,l,hmax);
        end
    end
    
end