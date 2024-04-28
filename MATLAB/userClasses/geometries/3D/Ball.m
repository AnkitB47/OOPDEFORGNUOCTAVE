classdef Ball < grid3D  
    methods(Access = public)
        function obj = Ball(varargin) 
            obj = obj@grid3D(); 
            obj.ball(varargin{:})
        end
    end
end