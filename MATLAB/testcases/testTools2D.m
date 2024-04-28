
% Do a 2D symmulation by using only low level (user) classes
% We test the typical workflow of solution of an 1D BVP using only classes
% Interval and LagrangeX1D. No pde object will be used.

DOM = UnitSquare(0.025);

FEM = Lagrange12D;
DOM.makeBoundaryMatrix(DOM.dirichletBC('-0.10'));

[K,~,F] = FEM.assema(DOM,1,0,10);
C = FEM.convection(DOM,[-10;-5]);
[~,~,H,R] = FEM.assemb(DOM);

s  = 100;
A = K+C+s*(H'*H);
b = F + s*H'*R;

u = A\b;

DOM.plot(u)




