classdef Lagrange21D < finiteElements1D
    % Class definition of P2 Elements for 1D space, here called
    % lagrange21D.
    % (C) 2015 Uwe Prüfert
    
    properties(Constant)
        % Define the elements stiffness matrices:
        % Stiff element
        S = reshape(1/3*[7 -8 1
                 -8 16 -8
                  1 -8 7],9,1);
        % Mass element
        M = reshape(1/30*[4 2 -1
                   2 16 2
                  -1 2 4],9,1);
              
        % Convection element      
        C = reshape([-0.5 0 0.5;...
                      0   0 0;...
                     -0.5 0 0.5],9,1);       
        % RHS element
        F = 1/6*[1 4 1]';  
        
        % Index of neighbours
        % We insert the additional point at number two and re-number the
        % left points to number three. Cf. the book of H.R.Schwarz. p.65.
        % This must be conform with the definition in the method grid1D.extendMesh
        % cf. obj.t =  ...
        idx  = 1:3;  
        assembleStep = 2;
    end
    
    methods(Static,Access = protected)
        % !!!! Dummies !!!!
        function fluxThrougElementEdges= fluxThroughEdges(gridObj,u,c)
        end
        function normFminusau = localErrorL2(gridObj,a,f) 
        end
        function jumps = fluxJumps(gridObj,fluxThroughElementEdges,order)
        end
    end 
    
    methods(Static,Access=public)  
        function [idxvec0,idxvec1,idxvec2] = makeIndex(idx,np)
            % Prepare the index sets for space call in abstract class.
            idxvec0 = reshape(idx,1,np*3);
            idxvec1 = reshape([idx;idx;idx],1,np*9);
            idxvec2 = reshape([idx(1,:);idx(1,:);idx(1,:);...
                                idx(2,:);idx(2,:);idx(2,:);...
                                idx(3,:);idx(3,:);idx(3,:)],1,9*np);     
        end
    end
    
    methods(Static,Access=protected)  
        function [J] = makeJ(gridObj) 
            scl = superclasses(gridObj);
            for k = 1:length(scl)
                isgrid = strcmp(scl{k},'gridd');
                if isgrid
                    break;
                end
            end
            if ~isgrid
                throw(obj.wrongClass);
            end
            %The first nElement entries are h =  x_i-x_i-1 from the global
            % mesh, the rest are entries from the extended mesh.
            J = gridObj.p(2:end)-gridObj.p(1:end-1);            
            J = J(1:gridObj.nElements);            
        end
    end 
end