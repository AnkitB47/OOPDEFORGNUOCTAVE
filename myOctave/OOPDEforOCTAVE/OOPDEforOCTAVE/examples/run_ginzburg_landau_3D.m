function run_ginzburg_landau_3D
ginzburglandau = GinzburgLandau();
ginzburglandau.grid = Ball(4*pi,0.25);
ginzburglandau.fem = Lagrange13D();
ginzburglandau.grid.plot
 
NT = 30;

u10 = @(x1,x2,x3) cos(x2);
u20 = @(x1,x2,x3) sin(x1);
 
x1 = ginzburglandau.grid.x;
x2 = ginzburglandau.grid.y;
x3 = ginzburglandau.grid.z;

ginzburglandau.y = [u10(x1,x2,x3)';u20(x1,x2,x3)'];

ginzburglandau.grid.plotFaces(u10(x1,x2,x3)');drawnow


ginzburglandau.time = linspace(0,NT,20*NT+1);
ginzburglandau.r =1;% 1.779077550000000; 
% ginzburglandau.r = 3; 

ginzburglandau.grid.makeBoundaryMatrix(ginzburglandau.grid.neumannBC('0'));
ginzburglandau.initialize
tic
ginzburglandau.solve
toc

%%
ginzburglandau.y = ginzburglandau.u1;% .^2-ginzburglandau.u2.^2;



 
%%
clf
ginzburglandau.animation('EdgeColor','none',...
    'Commands',...
    'colormap jet',...
    'view(3)',...
    'grid on')

end
 
