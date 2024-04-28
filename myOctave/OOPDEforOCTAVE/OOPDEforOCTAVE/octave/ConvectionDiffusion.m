classdef ConvectionDiffusion < pde & plotUtilsTimeDependent
       
    methods(Access = protected)
        function val = df(obj,~,y)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
           val = obj.A*y + obj.b;
        end
        
        function val = jacobian(obj,~,~)
            val = obj.A;
        end
    end
    
    methods(Access = public)        
         
        function initialize(obj,varargin)  
            initialize@pde(obj,varargin{:});            
            lambda = 1000000;           
            obj.A = -(obj.K+obj.M+obj.C+obj.Q+lambda*(obj.H'*obj.H));
            obj.b = lambda*obj.H'*obj.R+obj.G+obj.F;
        end
        
        function solve(obj,~)
            solve@pde(obj,'ODE15S'); 
        end
    end
end

