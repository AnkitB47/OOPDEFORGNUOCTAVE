function run_mountain  

    mountain = Mountain; 
   
    figure(1)
    clf
    mountain.plotFaces;hold on
    
    plot3(6.75,3.5,2.131,'.','MarkerSize',20)
   
    drawnow
    
    mountain.makeBoundaryMatrix(...
         mountain.neumannBC('0'),...
         mountain.neumannBC('0'),...
         mountain.neumannBC('0'),...
         mountain.neumannBC('0'),...
         mountain.neumannBC('0'),...    
         mountain.robinBC('1','1'));
 
 
 
    fem = Lagrange13D;

    pde = Elliptic;
    pde.grid = mountain;
    pde.fem = fem;

    f = source(mountain);
    figure(4)
    mountain.plotFaces(pde.grid.center2PointMatrix*source(pde.grid)')
    c = diffusity(mountain);

    pde.initialize(c,0,f)
    pde.solve
 
    figure(2)
    clf
    mountain.plotFaces(pde.y)
    colorbar
    drawnow

    
    % Postprocessing
    sub2d = Rectangle(0,10,0,4,0.125);

    emb1 = embeddedGridObject(mountain,sub2d);
    emb2 = embeddedGridObject(mountain,sub2d);
    emb3 = embeddedGridObject(mountain,sub2d);
    emb4 = embeddedGridObject(mountain,sub2d);
    emb5 = embeddedGridObject(mountain,sub2d);
    
    emb1.rotateSubmanifold(1,0,0,pi/2);
    emb2.rotateSubmanifold(1,0,0,pi/2);
    emb3.rotateSubmanifold(1,0,0,pi/2);
    emb4.rotateSubmanifold(1,0,0,pi/2);
    emb5.rotateSubmanifold(1,0,0,pi/2);
 
    emb1.moveSubmanifold(0,1,0)
    emb2.moveSubmanifold(0,3,0)
    emb3.moveSubmanifold(0,5,0)
    emb4.moveSubmanifold(0,7,0)
    emb5.moveSubmanifold(0,9,0)
    
    emb1.setData(pde.y)
    emb2.setData(pde.y)
    emb3.setData(pde.y)
    emb4.setData(pde.y)
    emb5.setData(pde.y)
   
    figure(3)
    clf
    emb1.plot('LineStyle','none')
    emb2.plot('LineStyle','none')
    emb3.plot('LineStyle','none')
    emb4.plot('LineStyle','none')
    emb5.plot('LineStyle','none')
    drawnow
%    
    
end

function val = source(grid)
    mdpts = grid.midpts;
    x = mdpts(1,:);
    y = mdpts(2,:);
    z = mdpts(3,:);
    val = zeros(1,grid.nElements);
    
    val(sqrt((x-5).^2+(y-5).^2+(z-1).^2)<.5) = 1e3;
end

function val = diffusity(grid)
    val = 3*ones(1,grid.nElements);
    val(grid.t(5,:)==2) = 20;
end
