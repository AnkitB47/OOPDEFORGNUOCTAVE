% Test script for solving the Schankenberg problem in 3D.
% For further remarks on the problem see the OOPDE Quickstart Guide.




% Define a bar. The mesh will have 16.000 nodes. 
kc=sqrt(sqrt(2)-1);
lx=4*pi/kc;
ly=4*pi/(sqrt(3)*kc);

% Call constructor
system = Schnakenberg();
system.grid = Bar(...
    linspace(-lx,lx,40),...
    linspace(-ly,ly,20),...
    linspace(-ly,ly,20)); 

system.fem =  Lagrange13D();

 

% Set lambda = 1;
system.lambda = 1.0;

% For compution the initial condition we need an "x".
x = system.grid.p(1,:); 

system.y = [system.lambda+cos(kc*x)'
            1/system.lambda-cos(kc*x)'];

% For time integration we use (0,1000) 
system.time = [0,1000];

% Call initialize w/o arguments.
system.initialize();

% Some ODE options
system.odeOptions = odeset('Stats','on');
system.solverOptions.solverTol = 1e-3;

% Use the steadyState solver. It iterates until ||df|| < solverTol. Be 
% patient. It can take 10 minutes.
system.solve('STEADYSTATE');

%
% % Solve the stationary problem.
% system.solverOptions.solverTol = 1e-6;
% system.solve('STATIONARY');
 
% Some visualization
figure(1) % Isofaces plot
clf;
system.grid.plotIso(system.u(:,end),[0.7  2  5]);
figure(2) % Plot faces ("Look from outside") 
clf;
system.grid.plotFaces(system.u(:,end),'LineStyle','none');
figure(3) % Cut slices ("Look inside")
clf;
system.grid.plotSlices(system.u(:,end),[-10 0 10],0,0,'LineStyle','none');



 
 



 