function run_burgers
    % Test program for burgers class.
    % It solves the viscous Burgers eqaution. Viscosity is the olny 
    % parameter for initialize. 
    % We use homogeneous Dirichlet boundary data on left and homogenious
    % Neumann boundary data on right boundary.
    % Initial value is classical piecewise linear hat function with y0(0) = 0,
    % y0(pi/2) = 1 and y0(pi) = 0. Try also y0 = sin(2x) by redefine the
    % function y0 at the end of the code.
    
    % Be aware that the spatial grid is refined properly.
    % After solving, a time-space plot and an animation will be shown.
    burgers = Burgers();
    burgers.fem = Lagrange11D;
    
    % Initialize spatial mesh. In 1D, the only possible domain is an interval.
    % Use Interval class to define the grid.    
    burgers.grid = Interval([0 pi]);
         
    % Refine spatial mesh seven times uniformly
    burgers.grid.refineUniformly(5);
    
    % Set up boundary condition

    burgers.setBoundaryConditions(...        
        'Dirichlet','0.0',...
        'Neumann','0');
    
    % Set y(0) = y0(x)
    burgers.y =  y0(burgers.grid.x);
    
    % Initialize burgers. Parameter is diffusity
    burgers.initialize(1e-2);
    
    % Chose time to solve, here [0,10] with 25 checkpoints. Try also 50,
    % 100, 200 checkpoint. For EULERI this determines also the time step
    % size as dt = 3/25, dt = 3/50 etc.
    burgers.time = linspace(0,.3,125);    
    
    % Solve with implicite Euler method by option 'EULERI'.
    % Try also 'BDF2', 'ODE15S' or 'ODE23S'.
    burgers.solve('BDF2');

    % Plotting in time-space style and animation over time.
    figure(1)
    burgers.animation('Color','k','Command','grid on')
    figure(2)
    burgers.plotTimeSpace('LineStyle','none','Command','axis equal');
    view(30,25)
end

% Initial value: sin(2x). 
function val = y0(x)
%     val = sin(x);
    val = 2*x/max(x);
    val(x>pi/2) = -2*(x(x>pi/2)-pi/2)/max(x)+1;
    val = val(:);
end
 
