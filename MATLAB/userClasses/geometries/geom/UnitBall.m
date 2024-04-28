classdef UnitBall < grid3D  
    methods(Access = public)
        function obj = UnitBall(varargin) 
            obj = obj@grid3D();         
            obj.unitBall(varargin{:})
        end
    end
end