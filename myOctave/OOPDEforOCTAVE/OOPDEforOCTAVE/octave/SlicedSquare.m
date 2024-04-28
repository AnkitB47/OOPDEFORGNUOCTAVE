classdef SlicedSquare < grid2D
    %slicedSquare Class definition
    %   This class defines a special geometry object.
    %   The domain is   [-1 1]^2 \{(-0.55 -0.45)x(-1,0) U (0.45 0.55)x(-1,0)} 
    %   i.e. a sliced square. 
    %   The class has a non-trivial constructor that uses distmesh2d
    %   functions to create the geometry object. In this sense, it is a
    %   demonstration how to us distmesh in connection with OOPDE. Note that
    %   we use her the adapted version of distmesh2d provided with OOPDE.
    
     
    
    methods(Access = public)
        function obj = SlicedSquare()
            fd=@(p) ddiff(ddiff(drectangle(p,-1,1,-1,1),...
                drectangle(p,-0.55,-0.45,-1,0)),...
                drectangle(p,0.45,0.55,-1,0)); 
            fixPoints1 = [linspace(-1,-0.55,15);-ones(1,15)]';
            fixPoints2 = [linspace(-0.45,0.45,20);-ones(1,20)]';
            fixPoints3 = [linspace(0.55,1,15);-ones(1,15)]';
            
            fh=@(p) 0.3+drectangle(p,-0.55,-0.45,-1,0)+drectangle(p,0.45,0.55,-1,0);
            [p,e,t] = distmesh2d(fd,fh,.05,...
                                [-1,-1;1,1],...
                                [-1,1;...
                                  1,1;...
                                  fixPoints1;...
                                  fixPoints2;...
                                  fixPoints3;...
                                -0.55,0;...
                                -0.45,0;...
                                0.45,0;...
                                0.55,0]);
            
            obj.p = p';
            obj.t = [t';ones(1,size(t',2))];
            obj.e = e';
            
           
            indx = find(obj.p(1,:)>-0.6&obj.p(1,:)<=-0.4&obj.p(2,:)>-1&obj.p(2,:)<=-0.2);          
            for k = 1:length(indx)
                [~,edge2] = find(obj.e(1:2,:)== indx(k));                 
                obj.e(5,edge2)=2;
            end
            
            indx = find(obj.p(1,:)>0.4&obj.p(1,:)<0.6&obj.p(2,:)>-1&obj.p(2,:)<=-0.2);          
            for k = 1:length(indx)
                [~,edge2] = find(obj.e(1:2,:)== indx(k));                 
                obj.e(5,edge2)=3;
            end
            
            indx = find(obj.p(1,:)>-0.6&obj.p(1,:)<-0.4&obj.p(2,:)>-0.2&obj.p(2,:)<=0.1);          
            for k = 1:length(indx)
                [~,edge2] = find(obj.e(1:2,:)== indx(k));                 
                obj.e(5,edge2)=4;
            end
            
            indx = find(obj.p(1,:)>0.4&obj.p(1,:)<0.6&obj.p(2,:)>-0.2&obj.p(2,:)<=0.1);          
            for k = 1:length(indx)
                [~,edge2] = find(obj.e(1:2,:)== indx(k));                 
                obj.e(5,edge2)=5;
            end            
        end
    end  
end

