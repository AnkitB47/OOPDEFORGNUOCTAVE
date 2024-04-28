%% A boundary value problem in $R^2$
%
%% First Example
%
% Define the pde. For an elliptic PDE
%
%
%
% $$\nabla \cdot (c(x) \nabla u) + a u = f$$
%


%%
% We use Elliptic class.
clc
elliptic = Elliptic();



%%
% Define an unit circle
clc
elliptic.grid = UnitCircle(0.125);


%%
% For linear elements on 2D domain use Lagrange12D
clc
elliptic.fem = Lagrange12D();


%%
%Refine the mesh two times uniformly
elliptic.grid.refineUniformly(2);

%%
% Define Dirichlet boundary conditions

elliptic.setBoundaryConditions(...
    'Dirichlet','sin(4*s).*sin(7*s).*sin(13*s)')
%%
% Call initialize. Parameters are diffusity c = 1, a = 0 and source
% f = 0.
c{1} = 1;

c{2} = ones(1,elliptic.grid.nElements);
c{3} = ones(1,elliptic.grid.nPoints);

c{4} = [1 0
        0 1];

c{5} = [ones(1,elliptic.grid.nPoints),ones(1,elliptic.grid.nPoints)];
c{6} = [ones(1,elliptic.grid.nPoints)
    ones(1,elliptic.grid.nPoints)];

c{7} = [ones(1,elliptic.grid.nElements),ones(1,elliptic.grid.nElements)];
c{8} = [ones(1,elliptic.grid.nElements)
    ones(1,elliptic.grid.nElements)];

c{9} = [1 0.3
    0.3 1];

c{10} = [ones(1,elliptic.grid.nElements) 0.3*ones(1,elliptic.grid.nElements)
     0.3*ones(1,elliptic.grid.nElements) ones(1,elliptic.grid.nElements) ];

c{11} = [ones(1,elliptic.grid.nPoints) 0.3*ones(1,elliptic.grid.nPoints)
     0.3*ones(1,elliptic.grid.nPoints) ones(1,elliptic.grid.nPoints)];

c{12} = {'1' '0.3'
         '0.3' '1'};

for k = 1:12

elliptic.initialize(c{k},0,0)

%%
% Solve linear proeblem. Use Algebraic multigrid solver by 'AMG'.
% Note that AMG needs   Ilupack. If Ilupack is not availiable,
% try 'LINEARGAUSS' or 'LINEAR' solver.

elliptic.solve('LINEARGAUSS')

%%
% Create (if there is no figure 1) and clear the figure.
% Plot the result. Since the mesh is rather fine, we use "LineStyle" =
% 'none' to supress printing the black edges of the triangles.
figure(1)   ;    clf      ;
elliptic.plot('LineStyle','-')
colormap cool
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
%% Second example.
% PDE on unit square with P2 elements and mixed Boundary
% conditions





%%
% This class uses stiff-spring technique to involve Dirichlet boundary
% conditions. For an example using Lagrange multiplier technique
% use
%

clc
elliptic = EllipticLagrange()

%
%%
% Define an L-shaped domain.
clc
elliptic.grid = UnitSquare(0.125)

%%
% For using  P2 elements you must first extend the mesh
%
% Note: When using linear elements, don't call extendMesh. Otherwise, an
% exeption will be thrown.
elliptic.grid.extendMesh;

% Use now P2 elements
%
clc
elliptic.fem = Lagrange22D()



%%
% The L-Shaped object created by Lshape class has six boundary
% segments. You can check how many boundary segments a geometry has
% by calling elliptic.grid.nBoundarySegments.
% The method elliptic.grid.identifyBoundarySegment provides a figure
% where the boundary segments can be identified by color  and number.

%%
% Three Dirichlet, one Neuman BC.

elliptic.setBoundaryConditions(...
    'Dirichlet','sin(pi*s).^3',...
    'Dirichlet','y',...
    'Dirichlet','x',...
    'Neumann','0')


%%
% Call initialize. Parameters are diffusity c = 0.1, a = 0 and source
% f = 10.
elliptic.initialize(1,0,10);

%%
% Solve linear proeblem. Use Algebraic multigrid solver by 'AMG'.
% Note that AMG needs   Ilupack. If Ilupack is not availiable,
% try 'LINEARGAUSS' or 'LINEAR' solver.
elliptic.solve('LINEARGAUSS');

%%
% Plot the result. Since the mesh is rather fine, we use "LineStyle" =
% 'none' to supress printing the black edges of the triangles.
figure(1)  ;    clf;
elliptic.plot('LineStyle','-');
colormap cool


%% Third Example Dirichlet BCs via Lagrange Multiplier
% Use Lagrange Multiplier to involve Dirichlet BCs
clc
elliptic = EllipticLagrange()

%
%%
% Define an L-shaped domain.
clc
elliptic.grid = Lshape(0.125)

%%


%%
% For linear elements on 2D domain use Lagrange12D
clc
elliptic.fem = Lagrange12D()

%%
% Note when using linear elements, don't call extendMesh. Otherwise, an
% exeption will be thrown.

%%
%Refine the mesh two times uniformly
elliptic.grid.refineUniformly(4);

%%
% The L-Shaped object created by Lshape class has six boundary
% segments. You can check how many boundary segments a geometry has
% by calling elliptic.grid.nBoundarySegments.
% The method elliptic.grid.identifyBoundarySegment provides a figure
% where the boundary segments can be identified by color  and number.

%%
% Since we want to have the same boundary condition on all six
% segments, we can use the simple call
elliptic.setBoundaryConditions(...
    'Dirichlet','1',...
    'Dirichlet','y',...
    'Neumann','0',...
    'Neumann','0',...
    'Dirichlet','x',...
    'Dirichlet','1'...
    );


%%
% Call initialize. Parameters are diffusity c = 0.1, a = 1 and source
% f = 10.
elliptic.initialize(1,0,10);

%%
% Solve linear proeblem. Use Algebraic multigrid solver by 'AMG'.
% Note that AMG needs   Ilupack. If Ilupack is not availiable,
% try 'LINEARGAUSS' or 'LINEAR' solver.
elliptic.solve('AMG');

%%
% Create (if there is no figure 1) and clear the figure.

figure(1)
%%
% Plot the result. Since the mesh is rather fine, we use "LineStyle" =
% 'none' to supress printing the black edges of the triangles.

elliptic.plot('LineStyle','-');
colormap cool



%% Fourth Example: Use OOPDE as Utility
%
% Use lowlevel methods to compute the matrices for the linear system
% and solve it explicitely by using bachslash.
% We don't use here a pde class object.

%%
% start with the domain, here a rectangle with rectangular elements.
clc
grid = RectangleR(0, 2*pi,0, pi)

%%
% Try to refine the mesh
grid.refineMesh;

%%
% Define Boundary Conditions: Here we declare only the type of boundary
% conditions for the four boundary segments.

grid.makeBoundaryMatrix(...
    grid.dirichletBC,...
    grid.dirichletBC,...
    grid.robinBC,...
    grid.robinBC)

%%
% Choose Bilinear elements (we are on a grid of rectangel type.)

fem = Bilinear2D

%%
% Assemblle matrices connected with the domain.

[K,M,F] = fem.assema(grid,1,1,0);

%%
% Determine index of boundary points for Boundary Segments No 1 & 2.

bp1 = grid.getBoundaryPointsIndexPerSegment(1);
bp2 = grid.getBoundaryPointsIndexPerSegment(2);

% compute the values at x anfd y coordimate

r = [sin(grid.x(bp1)),2*sin(-2*grid.y(bp2))] ;

%%
% Assemble boundary matrices, use data for Dirichlet, use simple constants
% for Robin BCs

[H,R] = fem.assembDataDirichlet(grid,1,r);
[Q,G] = fem.assembDataRobin(grid,1,0);

%%
% Soleve the linear system by using "\"
% Penalty is set to 10³
y = (K+M+1e3*(H'*H)+Q)\(F+1e3*(H'*R)+G);


%%
% Since there is no object derived from pde, we use lower lever plot
figure(1)   ;   clf;


grid.plot(y)

