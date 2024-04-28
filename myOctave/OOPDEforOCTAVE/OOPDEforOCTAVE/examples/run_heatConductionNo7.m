% Testprogram for heatConductionNo7 class
% Implements the problem from the Artificial Bee Colony (ABC-Algorithm) paper 

pde = heatConductionNo7();
pde.initialize();
pde.solve();
figure(1)
pde.plotTimeSpace();

colormap hot




