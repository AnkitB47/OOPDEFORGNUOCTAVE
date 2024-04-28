% Test for grid2DR and Bilinear2D 
% TODO Dirichlet BCs

grid = Rectangle(0, 2*pi,0, pi);
 
grid.refineMesh;
 
grid.makeBoundaryMatrix(...
    grid.robinBC,...   
    grid.robinBC,...
    grid.robinBC,...
    grid.robinBC);
  
fem = Lagrange12D;

[K,M,F] = fem.assema(grid,1,1,0);
 
 
bp1 = grid.getBoundaryPointsIndexPerSegment(1);
bp2 = grid.getBoundaryPointsIndexPerSegment(2);
bp3 = grid.getBoundaryPointsIndexPerSegment(3);
bp4 = grid.getBoundaryPointsIndexPerSegment(4);

bs1 = grid.getBoundaryElementsIndex(1);
bs2 = grid.getBoundaryElementsIndex(2);
bs3 = grid.getBoundaryElementsIndex(3);
bs4 = grid.getBoundaryElementsIndex(4);

xb = grid.point2CenterB(grid.x);
yb = grid.point2CenterB(grid.y);

g = [sin(xb(bs1)),...
    2*sin(2*yb(bs2)),...
    pi-0.5*xb(bs3),...
    yb(bs4)] ;


q = g;

% [H,R] = fem.assembDataDirichlet(grid,1,r);  

[Q,G] = fem.assembDataRobin(grid,100,100*g);


y = (K+M+Q)\(F+G);

 
grid.plot(y)