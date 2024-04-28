classdef Schnakenberg < pde
   % Class to solve the Schnakenberg problem. 
   % It assumes homogeneous Neumann BCs 
   % To control the solution use coefficents c, lambda and initial
   % conditions. 
    
    properties(Access = public)
        % The birfucation parameter lambda. To store the solutions components
        % we use u and v. Here, all must be from double
        % class, not function_handles etc.
        lambda double = 3 % Initialized. 
        c double = 60
    end
    
    properties(SetAccess = private)
        u double; 
        v double;
    end
    
    methods
        function val = get.u(obj)
                val = obj.y(1:obj.grid.nPoints,:);
        end
        function val = get.v(obj)
                val = obj.y(obj.grid.nPoints+1:end,:);
        end
    end
    
    methods(Access = protected)
        function val = df(obj,~,y)
            % The right hand side of parabolic pde, divided into
            % linear, nonlinear and constant part. 
            % obj.A is a matrix, obj.b a vector and obj.nonlinear a vector
            % valued function, here with parameter y
            % y is the vector containing the solutions of both PDE's.
            % This is a rather abstract definition, the definition of A, b
            % takkes place in initialize while nonlinear is a method
            val = obj.A*y+obj.nonlinear(y)+obj.b;
        end
        
        function J = jacobian(obj,~,y)
            n = obj.grid.nPoints;            
            % u is the first, v the second part
            ul = y(1:n);
            % We need a row vector of u
            ul = ul(:);
            vl = y(n+1:end);
            % v must be ca column vector
            vl = vl(:);
            J = obj.A+[2*obj.M*sparse(1:n,1:n,ul.*vl),...
                obj.M* sparse(1:n,1:n,ul.^2);...
                -2*obj.M*sparse(1:n,1:n,ul.*vl),...
                -obj.M*sparse(1:n,1:n,ul.^2)];
        end
        
        function [value,isterminal,direction] = eventfun(obj,t,y)
            % Overwrite the event function from pde.
            % Stops the integrations when max|df| < obj.termination
            value =  norm(obj.df(t,y),Inf)-abs(obj.solverOptions.solverTol); 
            isterminal = 1;   % stopps the integration
            direction = 0;   
        end   
    end
    
    methods(Access = public) 
        function initialize(obj)
            % Boundary conditions:  
            % We assume that all BC are homogeneous Neumann.
            
            % Assemble stiffness for c = 1, mass for a = -1 and source for
            % lambda (**)
            [K,M,L] = obj.fem.assema(obj.grid,1,1,obj.lambda);              
            n = sparse(obj.grid.nPoints,1);
            % Create a ultimate sparse to fill obj.A.
            N = sparse(obj.grid.nPoints,obj.grid.nPoints);
            
            % The linear part of the system, first line equation for u,
            % second line equation for v. The diffusion coefficient
            % ("lambda") for equation two is
            % written explicitely into the equation. This works only for
            % scalar parameters, not for "functions" c(x)
            obj.A = -[K+M     N
                       N   obj.c*K];
            
            % Constant part. 
            obj.b = [n
                     L]; % (**)
                 
            % Mass matrix to evaluate the nonlinearity     
            obj.M = obj.mass;    
            % Correct size of D matrix  
            obj.D = [obj.mass N
                N obj.mass];
            obj.initialized = true;
        end 
        
        function solve(obj,varargin)
            % Overwritten solve method.
            solve@pde(obj,'ODE15S');            
        end
        
        function val = r(obj,t,y)
            val = obj.df(t,y);
        end   
    end
    
    methods(Access = private)
        function val = nonlinear(obj,y)
            % Number of points to find the position of u and v in y
            n = obj.grid.nPoints;            
            % u is the first, v the second part
            ul = y(1:n);
            % We need a row vector of u 
            vl = y(n+1:end);             
            % Evaluate  u^2 v and -u^2 v (in weak form)   
            val =  [obj.M*(ul.^2.*vl)
                   -obj.M*(ul.^2.*vl)];            
        end
    end    
end

