classdef Room < grid3D
    methods(Access = public)
        
        function obj = Room(h)            
            %Room Class constructor for the room model with door and
            %window.
            %
            % Defines a 3D model of a romm with one (partially glassified)
            % door and a window. It has four boundary segments, modelling
            % different types of materials 
            % * 2 Concrete wall  inside building
            % * 3 Concrete walls street-side
            % * 4 Wood (door)
            % * 5 Glass (window/streetside) 
            % * 6 Glass parts of the door(inside building)
            % Segment 1 is the floor with the floor heating system.
            % 
            % Constructor calls
            % 
            %   g = Room            % w/o argument 
            %
            % or
            %
            %   g = Room(hmax)      % w hmax argument.
            % 
            % Default hmax = 0.25. If hmin is given as argument, h is set
            % by h = min(0.25,hmax).
            %
            
            switch nargin
                case 0
                    h = 0.25;
                case 1
                    h = min(h,0.25);
            end
            
            obj.bar(0, 3.6, 0, 4, 0, 2.6,h);

            % Redefine all boundary segments
            obj.e(end,:) = 1;
            
            % Define windows and door
            % First compute midpoints of boundary elements
            xm = obj.point2CenterB(obj.x);
            ym = obj.point2CenterB(obj.y);
            zm = obj.point2CenterB(obj.z);
            
            % Define walls inner
            obj.e(end,zm>0) = 2;
            % Define walls outer
            obj.e(end,zm>0&ym==4) = 3;
            % Door:
            obj.e(end,(1<xm&xm<=2.1&zm<2.1&ym==0)) = 4;
            % Window:
            obj.e(end,(1.15<xm&xm<=2.45&1<zm&zm<2.1&ym>0)) = 5; 
            % Window in Door
            obj.e(end,(1.2<xm&xm<=1.85&zm<1.85&0.5<zm&ym==0)) = 6; 
        end
        
        function plot(obj)
            %plot Method that plots the Room object.
            %
            % The object will be plotted nicely: Windows are transparent,
            % walls are gray, floor is black. 
            y = ones(1,obj.nEdges);
            clf;
            colormap gray;
            % colors for floor, door and make glass "transparent".
            y(obj.e(end,:) == 1) = 0;
            y(obj.e(end,:) == 4) = 0.7;
            y(obj.e(end,:) == 5) = NaN; 
            y(obj.e(end,:) == 6) = NaN;
            
            obj.plotFaces(y,'EdgeColor',[0.6 0.6 0.6]);
            
            % Nice borders of the room
            line([0 0],[0 0],[0 2.6],'Color',[0.4 0.4 0.4]);
            line([3.6 3.6],[0 0],[0 2.6],'Color',[0.4 0.4 0.4]);
            line([0 0],[4 4],[0 2.6],'Color',[0.4 0.4 0.4]);
            line([3.6 3.6],[4 4],[0 2.6],'Color',[0.4 0.4 0.4]);
            
            line([0 3.6],[0 0],[2.6 2.6],'Color',[0.4 0.4 0.4]);
            line([0 3.6],[4 4],[2.6 2.6],'Color',[0.4 0.4 0.4]);
            line([0 0],[0 4],[2.6 2.6],'Color',[0.4 0.4 0.4]);
            line([3.6 3.6],[0 4],[2.6 2.6],'Color',[0.4 0.4 0.4]);            
            
            colorbar off;
        end
    end            
end