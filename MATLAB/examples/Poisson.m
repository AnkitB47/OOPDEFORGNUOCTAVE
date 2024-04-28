classdef Poisson < pde
    methods(Access = protected)
        function  dy = df(obj,~,y)
                s = 1e3;
                dy = -(obj.K+...
                      s*(obj.H'*obj.H)+obj.Q)*y + ...
                      obj.F+s*(obj.H'*obj.R)+obj.G;
        end 
        function j = jacobian(obj,~,~)
            s = 1e3;
            j = -(obj.K+s*(obj.H'*obj.H)+obj.Q);
        end
    end  
end