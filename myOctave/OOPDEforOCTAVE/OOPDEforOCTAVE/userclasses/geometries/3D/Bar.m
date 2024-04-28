classdef Bar < grid3D    
     
    
    methods(Access = public)
        function obj = Bar(varargin)
            obj = obj@grid3D(); 
            obj.bar(varargin{:});
        end
    end
    
end