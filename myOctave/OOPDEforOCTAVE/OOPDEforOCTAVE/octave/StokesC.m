classdef StokesC < pde
    properties(Dependent = true, GetAccess = public, SetAccess = private)
        velocity_x
        velocity_y
        pressure
    end

    % getter methods
    methods
        function val = get.velocity_x(obj)
            val = obj.y(:,1:obj.grid.nPoints);
        end
        function val = get.velocity_y(obj)
            val = obj.y(:,1+obj.grid.nPoints:2*obj.grid.nPoints);
        end
        function val = get.pressure(obj)
            val = obj.y(:,1+2*obj.grid.nPoints:3*obj.grid.nPoints);
        end
    end

    methods(Access = public)
        function initialize(obj,nue)
            bcwall = obj.grid.dirichletBC(0);
            inflow = obj.grid.dirichletBC('10*sin(pi*s)');
            outflow = obj.grid.neumannBC(0);

            obj.grid.makeBoundaryMatrix(inflow,bcwall,outflow);
            [~,~,H,BInOut] = obj.fem.assemb(obj.grid);

            obj.grid.makeBoundaryMatrix(bcwall,bcwall,bcwall)
            [~,~,HW,Bwall] = obj.fem.assemb(obj.grid);

            [K,~,~] = obj.fem.assema(obj.grid,nue,0,0);

            N = sparse(obj.grid.nPoints,obj.grid.nPoints);

            Bx = - obj.fem.convection(obj.grid,[1
                                             0]);
            By = - obj.fem.convection(obj.grid,[0
                                             1]);

            C = - 0.25*obj.grid.hmax^2*K;

            s = obj.fem.stiffSpring(K+H'*H);

            obj.A = -[K+s*(H'*H), N, Bx';
                  N, K+s*(HW'*HW), By';
                  Bx, By, C];

            obj.b =  [s*H'*BInOut;
                      s*HW'*Bwall;
                      zeros(obj.grid.nPoints,1)];

            obj.initialized = true;
        end

        function plot(obj,varargin)
            clf
            xa=min(obj.grid.x); xb=max(obj.grid.x);
            ya=min(obj.grid.y); yb=max(obj.grid.y);
            subplot(2,2,1)
            quiver(obj.grid.x,obj.grid.y,obj.velocity_x(1,:),obj.velocity_y(1,:),3);
            axis([xa xb ya yb]);view(2)
            subplot(2,2,2)
            obj.grid.plot(obj.velocity_x(1,:),varargin{:});
            axis normal;view(2)
            subplot(2,2,3)
            obj.grid.plot(obj.velocity_y(1,:),varargin{:});
            axis normal;view(2)
            subplot(2,2,4)
            obj.grid.plot(obj.pressure(1,:),varargin{:});
            axis normal;view(2)
        end
    end

    methods(Access = protected)
        function dy = df(~,~,~)
            % A dummy
            dy = [];
        end
    end

end
