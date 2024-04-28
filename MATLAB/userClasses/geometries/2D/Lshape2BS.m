classdef  Lshape2BS < Lshape
    methods(Access = public)
        function obj = Lshape2BS(varargin)            
            obj@Lshape(varargin{:});
            if isempty(obj.e)
                return
            else
                indx3  = obj.e(5,:) == 3;
                indx4 = obj.e(5,:) == 4;
                obj.e(5,:) = 1;
                obj.e(5,indx3|indx4) = 2;	
            end
        end
    end
    
end
