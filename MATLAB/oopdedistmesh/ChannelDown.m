classdef ChannelDown < grid2D
    %ChannelDown Class definition
    %   This class defines a special geometry object.
    %   The domain is (0 3) x (-0.5 1) \{(0 1)x(-0.5,0) U (2 3)x(-0.5,0)}
    %   i.e. a sliced square.
    %   The class has a non-trivial constructor that uses distmesh2d
    %   functions to create the geometry object. In this sense, it is a
    %   demonstration how to us distmesh in connection with OOPDE. Note that
    %   we use her the adapted version of distmesh2d provided with OOPDE.

    methods(Access = public)
        function obj = ChannelDown()
            fd=@(p) ddiff(ddiff(drectangle(p,0,3,-0.5,1),...
                drectangle(p,0,1,-0.5,0)),...
                drectangle(p,2,3,-0.5,0));
            fixPoints1 = [linspace(0,1,20);-ones(1,20)]';
            fixPoints2 = [linspace(1,2,20);-ones(1,20)]';
            fixPoints3 = [linspace(2,3,20);-ones(1,20)]';


            fh=@(p) 0.3+drectangle(p,0,1,-0.5,0)+drectangle(p,2,3,-0.5,0);
            [p,e,t] = distmesh2d(fd,fh,.05,...
                                [0,-0.5;3,1],...
                                [0,1;...
                                  3,1;...
                                  fixPoints1;...
                                  fixPoints2;...
                                  fixPoints3;...
                                0,0;...
                                1,0;...
                                2,0;...
                                3,0]);

            obj.p = p';
            %obj.t = [t';ones(1,size(t',2))];
            obj.t = t';
            obj.e = e';

            % Marking boundary segments
            indx = find(obj.p(1,:)<(0+1e-8)&obj.p(2,:)>0&obj.p(2,:)<1); % adding tolerance 1e-8 as boundary segment performs over a range of values initially
            for k = 1:length(indx)                                       % rather than taking accurate value the first time, as the domain starts from[1,0)
                [~,edge2] = find(obj.e(1:2,:)== indx(k));                % so this will exclude the first x-cordinate of domain
                obj.e(5,edge2)=2;
            end

            indx = find(obj.p(1,:)>=3&obj.p(2,:)>0&obj.p(2,:)<1);
            for k = 1:length(indx)
                [~,edge2] = find(obj.e(1:2,:)== indx(k));
                obj.e(5,edge2)=3;
            end

        end

        % Function for giving the ChannelDowns Geometry parameters
        function [vertices] = getGeometryEdges(obj)
            %% Returns the outer geometry points of the sliced square before meshing.
            %% Input:
            %% Output:
            %%      vertices - size =
             vertices = [0 0;  1 0; 1 -0.5;
                                    2 -0.5;
                                    2 0;
                                    3 0;
                                    3 1;
                                    0 1; 0 0];
        end
    end
end

