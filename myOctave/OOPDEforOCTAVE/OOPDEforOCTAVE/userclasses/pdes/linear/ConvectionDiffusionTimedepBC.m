classdef ConvectionDiffusionTimedepBC < pde &  plotUtilsTimeDependent
    properties(Access=private)
        HH  
        MM  
    end
    
    properties(Access=public)        
        r
        f 
    end
    
    methods(Access = public)
        function  dy = df(obj,t,y) 
            dy = obj.A*y + obj.b +...
                obj.HH*obj.r(t,obj.grid.p(1,:))...
                + obj.MM*obj.f(t,obj.grid.p(1,:));
        end
        
        function J = jacobian(obj,~,~)            
            J = obj.A;
        end
    end
    
    methods(Access = public)
        function initialize(obj,d,c,b,a,f)            
            initialize@pde(obj,d,c,b,a,f);
            s = obj.fem.stiffSpring(obj.K+obj.M+obj.C);
            obj.A = -(obj.K+obj.M+obj.C+s*(obj.H'*obj.H)+obj.Q);           
            obj.b =  obj.G + obj.F; 
            obj.MM = obj.mass;
            obj.HH = s*obj.H'*obj.H;            
        end
    end  
end