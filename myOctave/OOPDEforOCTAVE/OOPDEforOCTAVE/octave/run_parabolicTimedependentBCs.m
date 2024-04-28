function run_parabolicTimedependentBCs
    % Call constructor
    pde = ConvectionDiffusionTimedepBC();
    
    % 1D domain
    pde.grid = Interval([0,2*pi],0.05);
    
    % P1 elements
    pde.fem = Lagrange11D(); 
    % or P2 elements
    %     pde.fem = lagrange21D();
    %     pde.grid.extendMesh()
    % Note: lagrange21D without extendMesh will lead to completely wrong
    % results.
    
    % Boundary conditions. Left homogeneous Neumann, right h u = 1.     
    pde.grid.makeBoundaryMatrix(...
                pde.grid.neumannBC('0'),...
                pde.grid.dirichletBC('1'));
    pde.initialize(1,1,0,0,0)

    % Set the boundary and source function 
    pde.r = @r;
    pde.f = @f; 
    pde.time = 0:0.5:30;
     
    x = pde.grid.p(1,:);
    pde.y = y0(x);
    
    
    % Call 'ODE15S' 
    pde.solve('BDF2')
    
    % Some visualization
    figure(1)
    cla;
    pde.plotTimeSpace();
     
    figure(2)
    cla;
    pde.animation();
end  

% Define y0 und r(t)
function y = r(t,x)
    y = 2*sin(t)*ones(size(x));
    y = y(:);
end
function y = f(t,x)
    y = 0*exp(-t)*sin(x);
    y = y(:);
end
function y = y0(x)
    y = zeros(size(x));
    y(x<=pi) = sin(x(x<=pi));
    y = y(:); 
end  
     
