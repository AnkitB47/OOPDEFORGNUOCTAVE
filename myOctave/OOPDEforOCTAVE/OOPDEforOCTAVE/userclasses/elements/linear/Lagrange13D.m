%% lagrange13D
%
% Class that implements Lagrange-1 elements in 3D.

%% Copy style
%
% handle

%% Inheritance
% lagrange13D < finiteElements3D < finiteElements < handle

classdef Lagrange13D< finiteElements3D
    
    %% Constant properties with Access = public
        %
        % * S1, ...,S9  (double) *vectors* that store the element stiffness matrices.
        % * C1,...,C3   (double) *vectors* that store the element gradient matrices
        % * M           (double) *vector* that stores the element mass matrix
        % * F           (double) vector that stores the source mass vector.
        % * idx         (double) vector that stores order of Points in element field.
        % * boundaryIndex (double) = [1:4]
    properties(Constant)
        %  stiff
        S1 =  reshape([1/6    -1/6     0     0
                      -1/6     1/6     0     0  
                       0     0     0     0
                       0     0     0     0] ,16,1);        
        
        S2 = reshape( [1/6     0    -1/6     0
                       0     0     0     0
                      -1/6     0     1/6     0
                       0     0     0     0],16,1);

        S3 = reshape( [1/6     0     0    -1/6
                       0     0     0     0
                       0     0     0     0
                      -1/6     0     0     1/6],16,1);
                  
                  
        % Note: Here we have to multiply the integral by the
        % factor of two, if we want so include the
        % symmetrie of intefgral here and not in the
        % assembling... cf.  method createMatrixEntries of
        % finiteElements3D class
        % But the integral is not symmetric as long as c is not symmetric.
        % So we do not multiply by 2.
      

        S4 =  reshape( [ 1     0    -1     0
                        -1     0     1     0
                         0     0     0     0
                         0     0     0     0 ] /6 ,16,1);
       
       
        S5 =  reshape( [ 1     0     0    -1
                        -1     0     0     1
                         0     0     0     0
                         0     0     0     0 ] / 6 ,16,1  );


        S6 =  reshape( [1     0     0    -1
                         0     0     0     0
                        -1     0     0     1
                         0     0     0     0] / 6,16,1);

                  
        % S7 to S9 represents the non symetric case. So we do not have to
        % multiply the integrals by 2 because we do not regard c as
        % symmetric
        
        
        S7 =  reshape( [1    -1     0     0
                         0     0     0     0
                        -1     1     0     0
                         0     0     0     0 ] / 6 ,16,1);
                     
        S8 =  reshape( [ 1    -1     0     0
                         0     0     0     0
                         0     0     0     0
                        -1     1     0     0 ] / 6 ,16,1  );


        S9 =  reshape( [1     0    -1     0
                         0     0     0     0
                         0     0     0     0
                        -1     0     1     0] / 6  ,16,1);
                  
 
       
       % Element Mass      
       
        M = reshape([1 / 60      1 / 120	 1 / 120	 1 / 120	
                     1 / 120	 1 / 60	     1 / 120	 1 / 120	
                     1 / 120	 1 / 120	 1 / 60	     1 / 120	
                     1 / 120	 1 / 120	 1 / 120	 1 / 60],16,1);  
 
        % convection 3x
        C1= reshape(1/24*[  -1 1 0 0
                            -1 1 0 0
                            -1 1 0 0
                            -1 1 0 0],16,1);
        
        
        C2 =  reshape(1/24*[ -1 0 1 0
                             -1 0 1 0
                             -1 0 1 0
                             -1 0 1 0],16,1);
 
         
        C3 =  reshape(1/24*[ -1 0 0 1
                             -1 0 0 1
                             -1 0 0 1
                             -1 0 0 1],16,1);
          
        % RHS vector  
        F = 1/24 * [1
                    1
                    1
                    1];
           
           
        % the points numbers, we have linear functions on
        % tetrahedrals, so we need   four points in this order:
        idx = 1:4;
    end
    
    properties(Constant = true, Access = protected)
        %% Constant properties with Access = protected
        %
        % * boundaryElements@finiteElements = lagrange12D
        %
        
        boundaryElements = Lagrange12D;
    end
    
    methods(Static,Access = public)   
        function [idxvec0,idxvec1,idxvec2] = makeIndex(idx,np)
            %% Overwritten method
            %
            % * makeIndex
            % * assemb
            %
            % Overwritten by dummies:
            %
            % * fluxThroughEdges
            % * ocalErrorL2
            % * fluxJumps
            
            
            % reshape index vectors for using sparse in assema
            idxvec0 = reshape(idx,1,np*4);
            idxvec1 = reshape([idx;idx;idx;idx],1,np*16);
            idxvec2 =  reshape([idx(1,:);idx(1,:);idx(1,:);idx(1,:);...
                idx(2,:); idx(2,:); idx(2,:);idx(2,:);...
                idx(3,:);idx(3,:);idx(3,:);idx(3,:);...
                idx(4,:);idx(4,:);idx(4,:);idx(4,:)],1,16*np);        
        end        
    end
    methods(Access=public,...
            Static = true)
        

        function [H,R] = assembDataDirichlet(gridObj,arg1,arg2)
            
            switch nargin
                case 1 % No data given: Assume homogenious Dirichlet Problem 
                    arg1 = 1;
                    arg2 = 0;
                case 2 % Case where only R is given. It implements h = r
                    arg2 = arg1;
                    arg1 = 1;
                case 3
                    % Nothig to do.
            end
            
            if isscalar(arg1) && isscalar(arg2)
                    h = arg1(ones(1,gridObj.nDirichletBoundaryPoints));                 
                    R = arg2(ones(gridObj.nDirichletBoundaryPoints,1));                    
                    H = sparse(1:gridObj.nDirichletBoundaryPoints,...
                        unique(gridObj.indexOfDirichletBoundaryPoints),...
                        h,...
                        gridObj.nDirichletBoundaryPoints,...
                        gridObj.nPoints);  
                    return
            else
                h = arg1;
                r = arg2;


                [bp,id,~] = unique(gridObj.indexOfDirichletBoundaryPoints);

                H = sparse(1:gridObj.nDirichletBoundaryPoints,...
                        bp,...
                        h(1:gridObj.nDirichletBoundaryPoints),...
                        gridObj.nDirichletBoundaryPoints,...
                        gridObj.nPoints);  
                R =  r(id)';
            end
        end
        function [Q,G] = assembDataRobin(gridObj,arg1,arg2)
            %% assembDataRobin 
            % IN:grid3Dpr,double,double OUT:double,double
            %
            % Assembles matrix Q and G from data vectors.
            % Arguments are the grid, the data g (Neumann BC) or q and g.
            % 
            % The data vectors must contain the value of q(x) and g(x)
            % respectivelly in the barycenters of the boundary elements,
            % triangle or sqares, respectivelly. 
            %
            % The order of the data must be the order of elements in obj.e.
            % of the Robin and/or Neumann baundary conditions.
            
            % If the boundary elements in the boundary  segments are not
            % increasing we resort the data.
            switch nargin
                case 2
                    g = arg1;  
                    q = zeros(1,gridObj.nRobinBoundaryElements);
                case 3
                    g = arg2; 
                    q = arg1;
 
                otherwise %
                    MException('Lagrange13D:WRONGNUMBERARGS',...
                        'Wrong number of arguments.').throwAsCaller;
            end
            if ~isa(gridObj,'grid3D')
                MException('Lagrange13D:WRONGARGS',...
                        'The first argument must be a grid3D object.').throwAsCaller;
            end
            if gridObj.nRobinBoundaryElements~=length(g)||gridObj.nRobinBoundaryElements~=length(q)
                gridObj.wrongFormat.throwAsCaller;
            end       
            
            %  indexOfTriagledFaces 
            indexOfTriagledFaces = gridObj.indexOfRobinBoundaryElements;
            
            % Compute the number of triangles in this group of boundary
            % elements. 
           
            numberOfTriangledFaces = length(indexOfTriagledFaces);
           
            % Problem: Faces are 2D embedded into the 3D space. 
            % Compute "J" by formula from elemetary geometry.
            
            % x-y-coordinates of the three face edge points...             
            p1 = gridObj.p(1:3, gridObj.e(1,indexOfTriagledFaces));
            p2 = gridObj.p(1:3, gridObj.e(2,indexOfTriagledFaces));
            p3 = gridObj.p(1:3, gridObj.e(3,indexOfTriagledFaces));

           
            % get element corner coordinates, define two vektors x y
            x1 = p2(1,:)-p1(1,:);
            x2 = p2(2,:)-p1(2,:); 
            x3 = p2(3,:)-p1(3,:);  
            
            y1 = p3(1,:)-p1(1,:);                
            y2 = p3(2,:)-p1(2,:);             
            y3 = p3(3,:)-p1(3,:);
            
            % Compute cross product X x Y
            z1 = x2.*y3 - x3.*y2;
            z2 = x3.*y1 - x1.*y3;
            z3 = x1.*y2 - x2.*y1;
            
            % compute ||z||
            J = sqrt(z1.^2+z2.^2+z3.^2);
           
            
            % Use elementary mass matrix from Lagrang12D class to compute
            % the 2 dimensional integrals
 
            
            Qe = reshape(Lagrange12D.M*(J.*q),1,9*numberOfTriangledFaces);
            ge = reshape(Lagrange12D.F*(J.*g),1,3*numberOfTriangledFaces); 

            
            % rearanging like in assema...
            indexOfTrianglePoints = gridObj.e(1:3,indexOfTriagledFaces);
 
            % 
            [indx0,indx1,indx2] = Lagrange12D.makeIndex(indexOfTrianglePoints,numberOfTriangledFaces) ;       

            Q = sparse(indx1,indx2,Qe,gridObj.nPoints,gridObj.nPoints); 
            G = sparse(indx0,1,ge,gridObj.nPoints,1); 
        end
    end
    
    methods(Access=public)        
        function [Q,G,H,R] = assemb(obj,gridObj)
            
            % Idee Daten aus gridObj.b gewinnen
            % dann assembDataRobin
            % bzw assembDataDirichlet berechnen.
            H = sparse(1,gridObj.nPoints);
            R = sparse(1,1);
            h = [];
            r = [];    
            
            Q = sparse(gridObj.nPoints,gridObj.nPoints);
            G = sparse(gridObj.nPoints,1);
            q = [];
            g = [];
            % If no Dirichlet BCs, do nothing
            
            
            for k = gridObj.indexOfDirichletBoundarySegments                
                bp = gridObj.getBoundaryPointsIndexPerSegment(k);
                p = [gridObj.x(bp);gridObj.y(bp);gridObj.z(bp)];
                [~,~,hval,rval] = gridObj.boundCoefficients(p,gridObj.b(:,k));  
                if isscalar(hval)
                    hval = hval(ones(1,size(bp,2)));
                end
                if isscalar(rval)
                    rval = rval(ones(1,size(bp,2)));
                end
                
                    
                h = [h,hval]; %#ok<AGROW>
                r = [r,rval]; %#ok<AGROW>
            end
           
            if ~(isempty(h)||isempty(r))
                [H,R] = obj.assembDataDirichlet(gridObj,h,r);
            end
                        
            for k = gridObj.indexOfRobinBoundarySegments 
                
                    mp = gridObj.midpointsOfBoundaryElements(k);                    
                    [qval,gval,~,~] = gridObj.boundCoefficients(mp,gridObj.b(:,k));
                    if isscalar(qval)
                        qval = qval(ones(1,size(mp,2)));
                    end
                    if isscalar(gval)
                        gval = gval(ones(1,size(mp,2)));
                    end
                
                q = [q qval];
                g = [g gval];
            end
            
            if ~(isempty(q)||isempty(g))
                [Q,G] = obj.assembDataRobin(gridObj,q,g);
            end      
        end 
        
        function [DX,DY,DZ] = gradientMatrices(obj,grid)
            %%
            % * gradientMatrices IN: gridd OUT: double,double
            % 
            % Method that computes matrices DX DY, such that 
            %
            % grad u = [(DX*u)' ,(DY*u)' (DZ*u)']                       
            %
            % at the center of all triangles. 
            % DX and DY are nElements x nPoints matrices. Note that this is
            % not an exact computation but an approximation with only
            % linear convergence. 
            
            % We want to use makeJ, so it is not static...
            
            % The indices of the first, second and third points in each
            % triangle. 
            idx1 = grid.t(1,:);
            idx2 = grid.t(2,:);
            idx3 = grid.t(3,:);
            idx4 = grid.t(4,:);
            
            % p1 are the  values of the "first" point, p2 the values of
            % the "second" point, etc. (2 x nElement matrices)
            p1 = grid.p(:,(idx1));
            p2 = grid.p(:,(idx2));
            p3 = grid.p(:,(idx3));
            p4 = grid.p(:,(idx4));
            
            % Compute Jacobi determinat  
            J = obj.makeJ(grid);
            
            % x_2-x_1 etc. in every element
            x21 = (p2(1,:)-p1(1,:));
            x31 = (p3(1,:)-p1(1,:));
            x41 = (p4(1,:)-p1(1,:));
            
            % y_2-y_1 etc. in every element
            y21 = (p2(2,:)-p1(2,:));
            y31 = (p3(2,:)-p1(2,:));
            y41 = (p4(2,:)-p1(2,:));
             
            % z_2-z_1 etc. in every element
            z21 = (p2(3,:)-p1(3,:));
            z31 = (p3(3,:)-p1(3,:));
            z41 = (p4(3,:)-p1(3,:));
            
            %     | x21 y21 z21 |
            % A = | x31 y31 z31 |
            %     | x41 y41 z41 |
            %
            %  -1    1   |y31*z41-z31*y41  z21*y41-y21*z41 y21*z31-z21*y31 |
            % A  = ----- |z31*x41-x31*z41  x21*z41+z21*x41 z21*x31-x21*z31 |
            %      det A |x31*y41-y31*x41  y21*x41+x21*y41 x21*y31-y21*x31 |
            %  
            %      | A B C |  
            %    = | D E F |
            %      | G H I |
            
            
            A = (y31.*z41-z31.*y41)./J;    
            B = (z21.*y41-y21.*z41)./J;
            C = (y21.*z31-z21.*y31)./J;
            D = (z31.*x41-x31.*z41)./J;
            E = (x21.*z41-z21.*x41)./J;
            F = (z21.*x31-x21.*z31)./J; %#ok<*PROPLC>
            G = (x31.*y41-y31.*x41)./J;
            H = (y21.*x41-x21.*y41)./J;
            I = (x21.*y31-y21.*x31)./J;
           
            DX = sparse([1:grid.nElements,1:grid.nElements,...
                1:grid.nElements,1:grid.nElements],...
                [idx1 idx2 idx3 idx4],[-(A+B+C) A B C] ,...
                grid.nElements,grid.nPoints);
            
            
            DY = sparse([1:grid.nElements,1:grid.nElements,...
                1:grid.nElements,1:grid.nElements],...
                [idx1 idx2 idx3 idx4],[-(D+E+F) D E F] ,...
                grid.nElements,grid.nPoints);
            
            DZ = sparse([1:grid.nElements,1:grid.nElements,...
                1:grid.nElements,1:grid.nElements],...
                [idx1 idx2 idx3 idx4],[-(G+H+I) G H I] ,...
                grid.nElements,grid.nPoints);
        end
    end
    
    
    
    methods( Static,Access = protected)
        % !!!! Dummies !!!!
        function fluxThrougElementEdges= fluxThroughEdges(gridObj,u,c)
        end
        function normFminusau = localErrorL2(gridObj,a,f) 
        end
        function jumps = fluxJumps(gridObj,fluxThroughElementEdges,order)
        end
    end 
end

