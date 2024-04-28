%% Schankenberg problem in $R^2$.
%


% Test script for solving the Schankenberg problem in 2D.
% Note that p2path uses a much smarter non-linear solver.
%%
% Call constructor
schnakenberg = Schnakenberg();
%%
% Define grid
schnakenberg.grid = Circle(22,0.5);
%%
% Use P1 elements
schnakenberg.fem = Lagrange12D();
%%
% If you want to try P2 elements, use this code
% system.fem = lagrange22D();
% system.grid.extendMesh();
%%
% Set lambda. Try , 1.0 2.0 3.0 and 3.1 
schnakenberg.lambda = 1.0;
schnakenberg.c = 60;
%%
% Call initialitze, no arguments needed.
schnakenberg.initialize()
%%
% For the initial condition we need an "x".
x = schnakenberg.grid.x;
%%
% Define initial condition
schnakenberg.y = [schnakenberg.lambda+cos(sqrt(sqrt(2)-1)*x)'   
            1/schnakenberg.lambda-cos(sqrt(sqrt(2)-1)*x)'];
        
schnakenberg.y = [(schnakenberg.grid.x+schnakenberg.lambda*schnakenberg.grid.y)'
    (schnakenberg.grid.x)'];        
%%
% Some ODE options
schnakenberg.odeOptions = odeset('Stats','on');
schnakenberg.solverOptions.solverTol = 1e-4;
%%
% Define time interval
schnakenberg.time = 0:500;
%
% Solve problem by time integration to obtain a initial guess for
% stationary solver.
schnakenberg.solve('ODE15S');

%%
% Animation
figure(1)
for k = 1:1:length(schnakenberg.time)    
   
    schnakenberg.grid.plot(schnakenberg.u(:,k),'LineStyle','none');
    view(2); 
    title(['t = ',num2str(schnakenberg.time(k))]);
    
    drawnow;
end


 