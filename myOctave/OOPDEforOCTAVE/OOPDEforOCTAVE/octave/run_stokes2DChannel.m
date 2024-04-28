
clc
stokes = StokesC();

stokes.grid = Channel(0.1);
stokes.grid.refineUniformly(1);
stokes.fem = Lagrange12D();

stokes.initialize(1)

stokes.solve('LINEARGAUSS')

fig = figure(1);

clf
stokes.plot('EdgeColor','none')
colormap hsv





