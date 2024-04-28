classdef heatConductionNo7 < pde & plotUtilsTimeDependent
    %heatConductionNo7 Class definition file.
    %   Class that implements parabolic PDE with time
    %   dependent Robin boundary condition
    
    
    properties(Constant = true)
        rhocp = 1000*2679
        lambda = 240
        d = 1
        T = 1000
        u0 = 980
    end
    
    properties(Access = public)
        a1 = 250
        a2 = 150
        a3 = 28
    end
    
    methods(Access = protected)
        function val =  df(obj,t,y)
            a = obj.boundary(t);
            obj.grid.makeBoundaryMatrix(...
                obj.grid.neumannBC('0'),...
                obj.grid.robinBC(num2str(a),num2str(a*298)));
            [obj.Q,obj.G,~,~] = obj.fem.assemb(obj.grid);           
            val = -(obj.A+obj.Q)*y+obj.G;
        end
        
        function val = jacobian(obj,~,~)
            val = -(obj.A+obj.Q);
        end
    end
    
    methods(Access = private)
        function val = boundary(obj,t)
            if t <= 90
                val = obj.a1;
            elseif t <= 250
                val = obj.a2;
            else
                val = obj.a3;
            end
        end
    end
    
    methods(Access = public)
        function initialize(obj,varargin)
            obj.fem = Lagrange11D;
            obj.grid = Interval([0,obj.d],0.01);
            obj.time = 0:1:obj.T;
            
            a = obj.boundary(0);
            obj.grid.makeBoundaryMatrix(...
                obj.grid.neumannBC('0'),...
                obj.grid.robinBC(num2str(a),num2str(a*298)));            
            initialize@pde(obj,obj.rhocp,obj.lambda,0,0,0);
           
            obj.A = obj.K;
            obj.y =  obj.u0*ones(obj.grid.nPoints,1);            
        end
        
        function solve(obj,varargin)
            solve@pde(obj,'ODE15S');
        end
    end
end

