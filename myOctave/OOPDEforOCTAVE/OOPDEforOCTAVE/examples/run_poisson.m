function run_poisson
    pde = Poisson();
    pde.fem = Lagrange12D();
    pde.grid = UnitSquare(0.1);
    pde.grid.refineMesh();
    pde.grid.makeBoundaryMatrix(pde.grid.dirichletBC('0'));
    pde.initialize(1,[0;0],0,10);
    pde.solve('LINEAR');
    pde.grid.plot(pde.y); 
end