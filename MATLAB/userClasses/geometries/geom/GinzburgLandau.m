classdef GinzburgLandau < pde & plotUtilsTimeDependent
    % Class that implements the Ginzburg-Landa equation
    
    properties(SetAccess = private)
        u1 double
        u2 double
    end
    
    properties(Access = public)
        r double;
    end
    
    methods(Access = protected)
        function val = df(obj,~,y)
            val = obj.A*y + obj.nonlin(y);  
        end
    end    
    
    methods(Access = public)
        function initialize(obj)
            [K,obj.M,~] = obj.fem.assema(obj.grid,1,1,0);
            [Q,~,H,~] = obj.fem.assemb(obj.grid); 
            l = 1e3;            
            N = sparse(obj.grid.nPoints,obj.grid.nPoints);
            obj.A = [-(K+Q+l*(H'*H))+obj.r*obj.M -obj.M
                      obj.M                      -(K+Q+l*(H'*H))+obj.r*obj.M];

            obj.D = [obj.M N
                     N obj.M];
            % To lazy to implement a jacobian, but give a pattern...    
            obj.pattern = obj.A~=0; 
            obj.initialized = true;      
        end
        
        function solve(obj)
            solve@pde(obj,'ODE15S')
            obj.u1 = obj.y(1:obj.grid.nPoints,:);
            obj.u2 = obj.y(1+obj.grid.nPoints:2*obj.grid.nPoints,:);
        end
    end
    
    methods(Access = private)
        function val = nonlin(obj,y)
            uu1 = y(1:obj.grid.nPoints);
            uu2 = y(1+obj.grid.nPoints:2*obj.grid.nPoints);
            
            val = -[obj.M*((uu1.^2+uu2.^2).*(-uu1-uu2))+obj.M*(((uu1.^2+uu2.^2).^2).*uu1)
                    obj.M*((uu1.^2+uu2.^2).*( uu1-uu2))+obj.M*(((uu1.^2+uu2.^2).^2).*uu2)];
        end
    end    
end

