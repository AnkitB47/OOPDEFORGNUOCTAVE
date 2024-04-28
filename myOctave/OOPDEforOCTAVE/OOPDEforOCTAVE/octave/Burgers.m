classdef Burgers < pde & plotUtilsTimeDependent
    % burgers
    % Defines a class for solving burgers equation. Superclass is pde.
    % (c) Uwe Pruefert 2013

     
    
    methods(Access=public)
        function initialize(obj,c)            
            initialize@pde(obj,1,c,1,0,0);            
            l = obj.fem.stiffSpring(obj.K+obj.C);
            obj.A = -(obj.K+l*(obj.H'*obj.H)+obj.Q);
            obj.b = l*obj.H'*obj.R+obj.G;
        end
    end

    methods(Access=public)
        function  dy = df(obj,~,y)
            dy = obj.A*y+obj.nonLin(y)+obj.b;
        end         
    end

    methods(Access=private)
        function y = nonLin(obj,y)
            y = -0.5*obj.C*(y.*y);
        end
    end
end