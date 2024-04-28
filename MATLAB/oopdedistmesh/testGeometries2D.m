 Testfile for all 2D geometries. Test error management later.

C = Circle(0.125);
C.refineMesh()
figure(1)
C.plot

T = DoubleT(0.25)
T.refineMesh;
figure(2)
T.plot

E = Ellipse(3,2,1,0.5,0.25);
E.refineMesh;
figure(3)
E.plot

H = HoleInPlane(0.125)
H.refineMesh;
figure(4)
H.plot

L = Lshape(0.125)
L.refineMesh;
figure(5)
L.plot

L2 = Lshape2BS(0.125)
L2.refineMesh;
figure(6)
L2.plot

T = Triangle(0.2)
T.refineMesh;
figure(7)
T.plot

O = Oval(0.125)
O.refineMesh;
figure(8)
O.plot

 %This will NOT work untill freeGeometry can handle intersections
##P = Part(0.125)
##P.refineMesh;
##figure(9)
##P.plot

 %Unit square, created with PolygonGeometry. No structured grid.
Poly = PolygonGeometry([0 0
                        1 0
                        1 1
                        0 1 ],0.2);
figure(10)
Poly.plot

SSq = SlicedSquare()
figure(11)
SSq.plot

##
Channel = PolygonGeometry([0 0;  1 0; 1 1.5;
                        0 1.5; 0 0],0.07);

Channel = PolygonGeometry([0 0
                           1 0
                           1 -0.5
                           2 -0.5
                           2 0
                           3 0
                           3 1
                           0 1],0.2)
figure(12)
Channel.plot


%ChannelDown
CD1 = ChannelDown()
figure
CD1.plot


us = UnitSquare(10);
figure;
us.plot();





