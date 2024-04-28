classdef grid3DR < gridd
    %grid2DR Class for 2D rectangle meshes
    %   Minimalistic implementation.
    
    properties(SetAccess = protected,...
            GetAccess = public)
        nPointsInElements double = 8;
    end
    
    properties(Access = protected)
        xgrid Interval
        ygrid Interval
        zgrid Interval
    end
    
    properties(Constant = true)
        spaceDimension double = 2;
    end
    
    methods(Access = public)
        function obj = bar(obj,x,y,z,h)
            
            switch nargin
                case 1
                    if isempty(obj.xgrid)||isempty(obj.ygrid)||isempty(obj.zgrid)
                        obj.xgrid = Interval(x);
                        obj.ygrid = Interval(y);
                        obj.zgrid = Interval(z);
                    end
                case 4
                    obj.xgrid = Interval(x);
                    obj.ygrid = Interval(y);
                    obj.zgrid = Interval(z);
                case 5
                     obj.xgrid = Interval(x,h);
                     obj.ygrid = Interval(y,h);
                     obj.zgrid = Interval(z,h);
            end
            nx = obj.xgrid.nPoints;
            ny = obj.ygrid.nPoints;
            nz = obj.ygrid.nPoints;
            
            x = obj.xgrid.x;            
            y = obj.ygrid.x;   
            z = obj.ygrid.x;  
            
            [x,y,z] = meshgrid(x,y,z);
            obj.p = [reshape(x,1,nx*ny*nz)
                reshape(y,1,nx*ny*nz)
                reshape(z,1,nx*ny*nz)];
            
             
             
            obj.t  = zeros(8,(nx-1)*(ny-1)*(nz-1));
            for k1 = 1:nz-1
                for k2 = 1:ny-1
                    for k3 = 1:nx-1
                    (k1-1)+(k2-1)+(k3-1)
                  [    1
                      ny+1
                      ny+2
                      2
                      nx*ny+1
                      nx*ny+ny+1
                      nx*ny+ny+2
                      nx*ny+2]
                      
                    
% % %                     (l-1)*(nz-1)+(k-1)*(ny-1)  
% % %                     1:ny-1
% % %                     (l-1)*(nz-1)+(k-1)*(ny-1)+(1:ny-1)
% % %                     obj.t(:,(l-1)*(nz-1)+(k-1)*(ny-1)+(1:ny-1)) = ...
% % %                         (l-1)*(nz)+(k-1)*(ny)+[ 1:ny-1
% % %                                                  ny+1:2*ny-1
% % %                                                  ny+2:2*ny
% % %                                                  2:ny
% % %                                                  1+nz*ny:ny*nz+ny-1
% % %                                                  ny*nz+nz+1:2*ny+nz*ny-1
% % %                                                  ny*ny+nz+2:2*ny+nz*ny
% % %                                                  2+nz*ny:ny*nz+nz  ]; 
                    end
                end
            end
            obj.t
        end
        
        function refineMesh(obj)
            obj.xgrid.refineMesh;
            obj.ygrid.refineMesh;
            obj.zgrid.refineMesh;
            obj.bar;
        end
      
        function [sidelength,area] = sideLengthAndArea(obj)
        end         
    end 
end