% Test for grid2DR and Bilinear2D 
% TODO Dirichlet BCs

grid = Part3D;
grid.refineMesh;
grid.plotFaces;

 
 
 
grid.makeBoundaryMatrix(...
    grid.robinBC(100,200),...   
    grid.robinBC(1000,0),...    
    grid.robinBC(100,200));
  
fem = Bilinear3D;

[K,M,F] = fem.assema(grid,1,0,0);
 
 
bp1 = grid.getBoundaryPointsIndexPerSegment(1);
bp2 = grid.getBoundaryPointsIndexPerSegment(2);
bp3 = grid.getBoundaryPointsIndexPerSegment(3);
 

bs1 = grid.getBoundaryElementsIndex(1);
bs2 = grid.getBoundaryElementsIndex(2);
bs3 = grid.getBoundaryElementsIndex(3);
 
xb = grid.point2CenterB(grid.x);
yb = grid.point2CenterB(grid.y);
g = [2*ones(size(xb(bs1))),...
    0*ones(size(xb(bs2))),...
    xb(bs3)] ;


q = g;

% [H,R] = fem.assembDataDirichlet(grid,1,r);  

[Q,G] = fem.assembDataRobin(grid,100*ones(size(g)),100*g);
% [Q,G,~,~] = fem.assemb(grid);

y = (K+M+Q)\(F+G);
clf
 
grid.plotFaces(y)

colorbar  