
% Do a 1D symmulation by using only low level (user) classes
% We test the typical workflow of solution of an 1D BVP using only classes
% Interval and LagrangeX1D. No pde object will be used.

DOM = Interval([-1 3],0.05);

FEM = Lagrange11D;
DOM.makeBoundaryMatrix(DOM.dirichletBC(0),DOM.dirichletBC('0'));

[K,~,F] = FEM.assema(DOM,1,0,1);
C = FEM.convection(DOM,'-1');
[~,~,H,R] = FEM.assemb(DOM);

s  = 100;
A = K+C+s*(H'*H);
b = F + s*H'*R;

u = A\b;

DOM.plot(u)




