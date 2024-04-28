%% lagrange22D 
%
% class that implements Lagrange-2 elements in 2D.
%
%% Copy style
% handle
% 
%% Inheritance 
% lagrange22D < finiteElements2D < finiteElements < handle
%


classdef Lagrange22D < finiteElements2D 
    % Written by Marcel Heim for his Diploma thesis. Some modification by
    % U.P. to include it "smooth" into the class lib OOPDE.
    properties(Constant,Access = public)
        
        %% Constant properties with Access = public
        %
        % * S1, S2 S3   (double) *vectors* that store the element stiffness matrices.
        % * C1 C2       (double) *vectors* that store the element gradient matrices
        % * M           (double) *vector* that stores the element mass matrix
        % * F           (double) vector that stores the source mass vector.
        % * idx         (double) vector that stores order of Points in element field.
        
        % definition of the stiffnes matrices
        S1 = reshape(1/6*[ 3  1 0 -4  0  0
                           1  3 0 -4  0  0
                           0  0 0  0  0  0
                          -4 -4 0  8  0  0
                           0  0 0  0  8 -8
                           0  0 0  0 -8  8
                  ],36,1);
              
            
              
        S2 = reshape(1/6*([ 6  1  1 -4  0 -4
                            1  0 -1 -4  4  0
                            1 -1  0  0  4 -4
                           -4 -4  0  8 -8  8
                            0  4  4 -8  8 -8
                           -4  0 -4  8 -8  8
            ]),36,1); 
              
        S3 = reshape(1/6*[ 3 0  1  0  0 -4
                           0 0  0  0  0  0
                           1 0  3  0  0 -4
                           0 0  0  8 -8  0
                           0 0  0 -8  8  0
                          -4 0 -4  0  0  8
                  ],36,1);
              
        % definition of the mass matrix      
        M = reshape(1/360*[  6 -1 -1  0 -4  0
                            -1  6 -1  0  0 -4
                            -1 -1  6 -4  0  0
                             0  0 -4 32 16 16
                            -4  0  0 16 32 16
                             0 -4  0 16 16 32
                  ],36,1);
        % definition of the convection matrices
                  
        C1 = reshape(...
              [-1/15           3/10           0             -7/30          -1/30           1/30    
               -3/10           1/15           0              7/30          -1/30           1/30    
                1/30          -1/30           0              0              1/15          -1/15    
                7/30          -7/30           0              0              2/15          -2/15    
                1/30           1/10           0             -2/15           4/15          -4/15    
               -1/10          -1/30           0              2/15           4/15          -4/15 
               ],36,1);  
           
     
        C2 = reshape(...
               [-1/15           0              3/10           1/30          -1/30          -7/30    
                 1/30           0             -1/30          -1/15           1/15           0       
                -3/10           0              1/15           1/30          -1/30           7/30    
                -1/10           0             -1/30          -4/15           4/15           2/15    
                 1/30           0              1/10          -4/15           4/15          -2/15    
                 7/30           0             -7/30          -2/15           2/15           0 
                 ],36,1);  
             
             
        % definition of the RHS element            
        F = 1/6*[0 0 0 1 1 1]';
     
        idx = 1:6; %changed. Order of Point in Points filed p. 
    end
        
    properties(Constant = true, Access = protected)
        %% Constant properties wirh Access = protected
        %
        % boundaryElements@finiteElements = lagrange21D
        %
        
        % For assembling the boundary condtion matrices, we need the class
        % of associated FE in 1D. Here we use lagrange21D.
        boundaryElements = Lagrange21D;     
    end
    
    properties(Constant,Access = public)       
        %% Constant properties with Access = public
        % * boundaryIndex@double = [1 6 2]
        boundaryIndex = [1 6 2]; 
        % Order of Points in boundary edge field. 
        % Here we have an extendet mesh, where additional points are added
        % in the last row. Note that the boundary element is lagrange21D.
    end
    
    methods(Static = true,Access=public)
        function [idxvec0,idxvec1,idxvec2] = makeIndex(idx,np)
            % Prepare the index sets for space call in abstract class.
            idxvec0 = reshape(idx,1,np*6);
            idxvec1 = reshape([idx;idx;idx;idx;idx;idx],1,np*6*6);
            idxvec2 =  reshape([idx(1,:);idx(1,:);idx(1,:);idx(1,:);idx(1,:);idx(1,:);...
                                idx(2,:);idx(2,:);idx(2,:);idx(2,:);idx(2,:);idx(2,:);...
                                idx(3,:);idx(3,:);idx(3,:);idx(3,:);idx(3,:);idx(3,:);...
                                idx(4,:);idx(4,:);idx(4,:);idx(4,:);idx(4,:);idx(4,:);...
                                idx(5,:);idx(5,:);idx(5,:);idx(5,:);idx(5,:);idx(5,:);...
                                idx(6,:);idx(6,:);idx(6,:);idx(6,:);idx(6,:);idx(6,:);]...
                                                  ,1,6*6*np);     
        end  
        
    end
    
    methods(Access=public)
        function [K,M,F] = assema(obj,gridObj,c,a,f)
            if isa(f,'double')&&length(f)==gridObj.nPoints
                [~,M,~] = assema@finiteElements2D(obj,gridObj,0,1,0);
                F = M*f(:);
                [K,M,~] = assema@finiteElements2D(obj,gridObj,c,a,0);                 
            else
                [K,M,F] = assema@finiteElements2D(obj,gridObj,c,a,f);
            end
        end
        
        function [DX,DY] = gradientMatrices(obj,grid)
            %%
            % * gradientMatrices IN: gridd OUT: double,double
            % 
            % Method that computes matrices DX DY, such that 
            %
            % grad u = [(DX*u)' ,(DY*u)']                       
            %
            % at the center of all triangles. 
            % DX and DY are nElements x nPoints matrices. Note that this is
            % not an exact computation but an approximation with only
            % linear convergence. 
            
            % We want to use makeJ, so it is not static...
            
            % The indices of the first, second and third points in each
            % triangle. Use "inner triangle" of extendet mesh
            idx1 = grid.t(4,:);
            idx2 = grid.t(5,:);
            idx3 = grid.t(6,:);
            
            % p1 are the  values of the "first" point, p2 the values of
            % the "second" point p3 the values of the third point 
            % of all triangles. (2 x nElement matrices)
            p1 = grid.p(:,(idx1));
            p2 = grid.p(:,(idx2));
            p3 = grid.p(:,(idx3));
            
            % Compute Jacobi determinat  
            J = obj.makeJ(grid);
            
            p12 = (p1(2,:)-p2(2,:))./J;
            p31 = (p3(2,:)-p1(2,:))./J;
            p23 = (p2(2,:)-p3(2,:))./J;
           
            DX = 4*sparse([1:grid.nElements,1:grid.nElements,1:grid.nElements],...
                [idx1 idx2 idx3],[p23 p31 p12] ,...
                grid.nElements,grid.nPoints);
            
            p21 = (p2(1,:)-p1(1,:))./J;
            p13 = (p1(1,:)-p3(1,:))./J;
            p32 = (p3(1,:)-p2(1,:))./J;
            
            DY = 4*sparse([1:grid.nElements,1:grid.nElements,1:grid.nElements],...
                [idx1 idx2 idx3],[p32 p13 p21] ,...
                grid.nElements,grid.nPoints);
        end  
    end
end
