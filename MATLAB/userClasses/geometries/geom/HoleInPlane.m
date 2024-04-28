classdef HoleInPlane < grid2D
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods
        function obj = HoleInPlane(hmax)
            % needs distmesh package
            if nargin == 0
                % default constructor
            else
            fd = @(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.4));
                pfix = [-1,-1;-1,1;1,-1;1,1]; 
                fh = @(p) min(sqrt( p(:,1).^2 + p(:,2).^2 ) , 1 );
                [p,e,t] = distmesh2d(fd,fh,hmax,[-1,-1;1,1],pfix);

                obj.p = p';
                obj.t = [t';ones(1,size(t',2))];
                obj.e = e';
            end
        end
    end
    
end

