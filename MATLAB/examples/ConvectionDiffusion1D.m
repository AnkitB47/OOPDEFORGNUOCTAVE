%% ConvectionDiffusion1D
%
%%
% This example shows how we define a specialized class that impelemnts only 
% one application.  
%
% We solve the pde (in weak) formulation
% 
% $$(d\dot u,v) +(c u_x,v_x) +(\vec b u_x,v)  + (a u,v) = (f,v) $$
%
% where  $u(0) = 1$ and   $u_x(1)=0$.
% The domain $\Omega$ is a subset of $R$. 
%
% This restiction allows us to set the grid and fem properties by the
% calss constuctor. 


%% Definition
% We derive the class from abstract classes pde and plotUtilsTimeDependent

classdef ConvectionDiffusion1D < pde & plotUtilsTimeDependent

% properties(Access = public)
%     fem
%     grid
% end
%% Override df and jacobian 
% The problems is a time dependent *linear* PDE. We can define df simply 
% (and rather formaly) by using the properties A and b, which we must
% define later. 
%
% Since the jacobian matrix of the problem is known, we can override
% jacobian method also by - here a very simple - user code.
    methods(Access = protected)
        function val = df(obj,~,y)            
           val = obj.A*y + obj.b;
        end
        
        function val = jacobian(obj,~,~)
            val = obj.A;
        end
    end
    
%% Public methods
%
% We uverride here the default constructor, the initialize and the solve
% method.
    
    methods(Access = public)
    
 %% Constructor
 % 
 % We override the default constructor by our owen code.
 % Since the domain cannot be any other thing than an interval, we
 % initialize it here. 
 %
 % Also the finite elements method is defined here and for ever by
 % Lagrange11D.
 %
 % Also the boundary conditions are defined here (and forever.) 
    
        function obj = ConvectionDiffusion1D
            obj.fem = Lagrange11D;
            
            obj.grid = Interval([0,1],0.001);          
            obj.setBoundaryConditions(...
               'Dirichlet', '1',...
               'Neumann','0');
        end
 %% Initialize
 % 
 % Intialize uss in its first line a call of the inherited intialize
 % method, giving all arguments to it. initialize@pde will compute all
 % matrices and  vectors needed for solving our problem.
 % lambda is a real number to mimicry a "stiff spring" to approximate the
 % Dirichlet boudary condition by a Robin boundary condition.
 % Last steps are the definition of matrix A and vector b.
        function initialize(obj,varargin)  
            initialize@pde(obj,varargin{:});            
            lambda = obj.fem.stiffSpring(obj.K+obj.M);            
            obj.A = -(obj.K+obj.M+obj.C+obj.Q+lambda*(obj.H'*obj.H));
            obj.b = lambda*obj.H'*obj.R+obj.G+obj.F;
        end
 %% solve
 %
 % To be sure to use the  right  solver
 % we override solve. Now, ODE15S  is the only
 % solver that can be used. 
        function solve(obj,~)
            solve@pde(obj,'ODE15S'); 
        end                
    end
    
%% check methods
% 
% For protecting the class for abuse, we restrickt to set fem only to
% Lagrange11D and grid to Interval objects.
% Its now not possible to use this class for 2D problems.
    methods(Static = true)
        function b = checkfem(fem)
            b = isa(fem,'Lagrange11D');
        end
        function b = checkgrid(grid)
            b = isa(grid,'Interval');
        end
    end
    
     
end

