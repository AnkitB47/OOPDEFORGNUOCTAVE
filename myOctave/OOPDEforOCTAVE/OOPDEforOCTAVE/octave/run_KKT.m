function run_KKT
% This is an example for solving 2D OCP. The desired state ist the solution
% of a Poisson Eqn with Dirichlet BCs and the Dirac source at 0 on the unit
% circle. The OCP has a Poissin Eqn with Neumann BCS!
    % construct a desired function yd as solution of a pde problem
    pde = Elliptic();  
    pde.grid = UnitCircle(0.0625);       
    % P1 elememts. Note that by the construction of the source for yd this
    % owrks NOT with P2 elements.
    pde.fem = Lagrange12D();   
    
    % create the comain: Unitcircle, refine twice and refine 5 time at
    % (0,0)
        
    for k = 1:5
        pde.grid.refineMesh(pde.grid.point2ElementIndex([0 0]));          
    end 
    
    % Poisson Eqn with Dirichlet BC
    pde.grid.makeBoundaryMatrix(pde.grid.dirichletBC('0')); 
    F = f(pde);
    pde.initialize(1,0,F);
    
    % solve the pde 
    pde.solve('LINEARGAUSS');
    subplot(2,2,4)   
    pde.grid.plot(pde.y);
    title('Desired state y_d')
    view(0,90)
    

    % define the kkt system, use copies of pde
    kktSystem = KKT();
    kktSystem.grid = pde.grid;
    kktSystem.fem = pde.fem;
    
    % change the boundary condition
    kktSystem.bcState = kktSystem.grid.neumannBC('0.0');    
        
    %  regularization lambda
    lambda = 1e-6;
    kktSystem.initialize(lambda,pde.y);  
    
    % solve KKT System
    kktSystem.solve('LINEARGAUSS'); 
    
    
    % Postprossesing
    figure(1)
    clf
    kktSystem.plot
   
    
    subplot(2,2,4)   
    pde.grid.plot(pde.y);
    title('Desired state y_d')
    view(0,90)
    
    
    % Model a point source with int_Omega |u| dx = 1;
    function y = f(p)
        y = zeros(size(p.grid.p(1,:)));
        ndx = p.grid.point2ElementIndex([0,0]);
        indx = unique(p.grid.t(1:3,ndx));
        [~,indxindx] = min(sum((p.grid.p(:,indx)-0).^2));
        y(indx(indxindx)) = 1e6;    
        y = y/(ones(size(y))*p.mass*abs(y)');        
    end
end