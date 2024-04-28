
clc
stokes.grid = ChannelDown();
stokes.grid.refineUniformly(2);
stokes.fem =  Lagrange12D();


stokes.initialize(1)

stokes.solve('LINEARGAUSS')

fig = figure(1);
set(fig,'Position',[400,200,800,600]);
clf
stokes.plot('EdgeColor','none')
colormap hsv





