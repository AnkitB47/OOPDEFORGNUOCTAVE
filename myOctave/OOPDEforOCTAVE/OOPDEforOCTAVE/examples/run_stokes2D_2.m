 
clc
stokes = Stokes2();
 
stokes.grid = channelWithCavity;  

stokes.fem = Lagrange12D();
 
stokes.initialize(1)

stokes.solve('LINEARGAUSS')  

fig = figure(1);
fig.Position = [400,200,800,600];
clf      
stokes.plot('EdgeColor','none')
colormap hsv
     
 


 
