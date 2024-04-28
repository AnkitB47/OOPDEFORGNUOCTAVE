%% lagrange03D
%
% Class that implements Lagrange-0 elements. The only method is mass.

%% Inheritance
% Lagrange03D < finiteElements3D

classdef Lagrange03D < finiteElements3D
    %class stub that implements L0 elements     
    
    properties(Constant)       
        idx@double = [];
    end
     
    
    methods(Static,Access = public)
        %% Static methods with Access = public
        %
        function M = mass(grid)
            %%
            % * mass IN:grid3D OUT:double
            %
            % Computes the mass matrix for Lagrange-0 elements. 
            %
            % Call:
            %
            %       M = Lagrange03D.mass(grid3D)
            %
            j =  Lagrange03D.makeJ(grid);            
            M = 1/6*sparse(1:grid.nElements,...
                1:grid.nElements,j);  
        end
    end
    
    methods( Static,Access = public)      
        % !!!! Dummies !!!!
        function makeIndex()
        end
    end
    methods( Static,Access = protected)
        function fluxThroughEdges()
        end
        function localErrorL2() 
        end
        function fluxJumps()
        end
    end
end

