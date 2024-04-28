classdef Channel < grid2D
    methods(Access = public)
        function obj = Channel(h)
            if nargin == 0
                h = Inf;
             end
             obj.polygonGeometry([0 0
                            1 0
                            1 -1
                            2 -1
                            2 0
                            4 0
                            4 1
                            0 1],h) ;

             % Create boundary segments for inflow an doutlet and no slip.

            obj.e(end,obj.e(end,:)==1) = 2;
            obj.e(end,obj.e(end,:)==2) = 2;
            obj.e(end,obj.e(end,:)==3) = 2;
            obj.e(end,obj.e(end,:)==4) = 2;
            obj.e(end,obj.e(end,:)==5) = 2;
            obj.e(end,obj.e(end,:)==7) = 2;

            obj.e(end,obj.e(end,:)==6) = 3;
            obj.e(end,obj.e(end,:)==8) = 1;
        end
    end
end


