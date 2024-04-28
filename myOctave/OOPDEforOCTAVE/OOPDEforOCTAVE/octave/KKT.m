classdef KKT < pde  
    properties(Access = public)
        % to allow the user to change the boundary condition 
        % of the state equation we make it public.
        bcState
    end
    properties(SetAccess = private)
        % propertise to store the state, the control and the adjoint
        % equation. Only methods from kkt should change the values.
        state
        control
        adjoint
    end
    
    methods (Access = protected)
        % These methods must be declared as protected. 
        function dy = df(obj,~,y)
            % More formal definition of the lienar pde system. A and be
            % must be implemented in initialize
            % The interface must be obj, time (not used here), and solution
            dy = obj.A*y+obj.b;            
        end  
        function J = jacobian(obj,~,~)
            % Since df is linear we know the Jacobian. it ist simply obj.A
            % The interface must be obj, time (not used here, and solution)
            J = obj.A;
        end
    end
    
    methods(Access = public)
        function initialize(obj,lambda,y_d)  
            % The most importand method of all
            % We know that the boundary condition for the adjoint is always
            % homogeneous. Hence, we have only non trivial boundary 
            % conditions for the state
                      
            obj.grid.makeBoundaryMatrix(obj.bcState);
            [Q,G,H,Ry] = obj.fem.assemb(obj.grid);
            % We have no source for the state equation but a source for the
            % adjoint. hence we can use the matrices for both equations and
            % F only for the adjoint.
            [K,M,F] = obj.fem.assema(obj.grid,1,0,...
                -obj.grid.point2Center(y_d));
            % The stiff-spring approximation of Dirichlet boundary
            % conditions (only needed if applyed
            l = obj.fem.stiffSpring(K+M);
            % Some sparse matrices/ vectors to fill the holes in the 
            % big matrix
            N = sparse(obj.grid.nPoints,obj.grid.nPoints);
            n = sparse(obj.grid.nPoints,1);
             
            % The pde system matrix. first line adjoint, second gradient
            % third state equation. Solution vector sate, control, adjoint
            obj.A = -[-obj.mass          N                K+M+l*(H'*H)+Q;... 
                      N                -lambda*obj.mass  -obj.mass;...
                      K+M+l*(H'*H)+Q   -obj.mass         N];              
            % The big rhs, F contains y-yd but no boundary vector
            % The state eqn. RHS has no source but both of the boundary 
            % vectors 
            obj.b = [F;...
                     n;...
                     l*H'*Ry+G];
            % To initialize y we set it to zero. This is very importand
            % since the solvers need the information that the number of
            % grid points in the mesh is nit equal with the number of dof.
            obj.y = zeros(obj.grid.nPoints*3,1);
            % At last, we set initialize true. Otherwise the solve method
            % returns an error.
            obj.initialized = true;
        end
        
        function solve(obj,varargin)
            % First call the solve from pde with options and then call
            % splitSolution method from kkt class
            if ~obj.initialized
                obj.notInitialized.throwAsCaller;
            end
            solve@pde(obj,varargin{:});
            obj.splitSolution
        end
            
        
        function plot(obj)  
            % plot works without splitSolution because we want to use a
            % loop.
            titles = {'State y' 'Control u ' 'Adjoint p'};
            for k = 1:3
                subplot(2,2,k)
                obj.grid.plot(obj.y((k-1)*obj.grid.nPoints+1:...
                    k*obj.grid.nPoints,end));
                view(2)
                title(titles{k});
            end
        end
    end
    
    methods(Access = private)
        function splitSolution(obj)
            % We split the solutions by simple indexing obj.y and using
            % obj.grid.nPoints
            obj.state = obj.y(1:obj.grid.nPoints,end);
            obj.control = obj.y(obj.grid.nPoints+1:2*obj.grid.nPoints,end);
            obj.adjoint = obj.y(2*obj.grid.nPoints+1:3*obj.grid.nPoints,end);
        end
    end    
end

