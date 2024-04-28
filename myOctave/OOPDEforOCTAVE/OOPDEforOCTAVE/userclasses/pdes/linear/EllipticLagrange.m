classdef EllipticLagrange < Elliptic     
    methods(Access=public)
        function initialize(obj,c,a,f)             
            initialize@Elliptic(obj,c,a,f);
           
            N = sparse(size(obj.H,1),size(obj.H,1));
          
            obj.A = -[(obj.K+obj.M)  obj.H'
                     obj.H           N];
                 
            obj.b = [obj.F + obj.G
                     obj.R];   
        end
        
        function solve(obj,varargin)
            % solve method, overwrites pde.solve. We must handle the
            % multipier. The matrix D must be overwriten.
            
            if isempty(obj.y)
                obj.y = [zeros(obj.grid.nPoints,1)
                    zeros(size(obj.R))];           
            end
            
            try
                solve@Elliptic(obj,varargin{:}); 
            catch ME
                ME.throwAsCaller;
            end
            obj.y = obj.y(1:obj.grid.nPoints,:);            
        end
    end
end