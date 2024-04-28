classdef DoubleT < grid2D
    %DoubleT Implements a grid for a double T geometrie.




    methods

        function obj = DoubleT(a1,a2,a3)
            %doubleT
            % Creates a grid for a double-t shaped domain
            % doubleT()
            % doubleT(scale)
            % Called with no argument, the base hal lenght one and the
            % height is three. One can scale the whole domain, but not
            % the proportions.
            % There are 12 boundary segments.
            % To identify boundary segment number call
            % obj.identifyBoundarySegment
            % (c) 2013 by Uwe Prüfert.

            switch nargin
                case 0
                    % default constructor
                     return
                case 1
                    % scale factor
                    hmax = a1;
                    scalex = 1;
                    scaley = 1;
                    if ~(isscalar(hmax)&&isa(hmax,'double'))
                        obj.wrongClass.throwAsCaller();
                    end
                case 2
                    if a1>0 && a2>0
                        scalex = a1;
                        scaley = a2;
                    else
                        ME = MException('GRIDD:POSITIVE',...
                             'Scale factors must be positive.');
                        ME.throwAsCaller();
                    end
                    hmax = Inf;
                case 3
                    if a1>0 && a2>0
                        scalex = a1;
                        scaley = a2;
                    else
                        ME = MException('GRIDD:POSITIVE',...
                             'Scale factors must be positive.');
                        ME.throwAsCaller();
                    end
                    if a3 > 0
                        hmax = a3;
                    else
                        ME = MException('GRIDD:POSITIVE',...
                             'hmax must be positive.');
                        ME.throwAsCaller();
                    end
                otherwise
                    obj.wrongNumberInputs.throwAsCaller();
            end

            obj.p = [0 1 2 3 0 1 2 3 1 2 1 2 0 1 2 3 0 1 2 3;...
                     0 0 0 0 1 1 1 1 2 2 3 3 4 4 4 4 5 5 5 5];

            obj.e = [1   2   3   4 8 7   10  12  15 16  20  19  18  17 13 14  11   9   6  5;...
                     2   3   4   8 7 10  12  15  16 20  19  18  17  13 14 11  9    6   5  1;...
                     0   1/3 2/3 0 0 0   1/3 2/3 0  0   0   1/3 2/3 0  0  0   1/3  2/3 0  0;...
                     1/3 2/3 1   1 1 1/3 2/3 1   1  1   1/3 2/3 1   1  1  1/3 2/3  1   1  1;...
                     1   1   1   2 3 4   4   4   5  6   7   7   7   8  9  10  10   10  11 12];

            obj.t = [1 2 2 3 3 4 6 7  9  10 11 12 13 14 14 15 15 16;...
                     2 6 3 7 4 8 7 10 10 12 12 15 14 18 15 19 16 20;...
                     5 5 6 6 7 7 9 9  11 11 14 14 17 17 18 18 19 19;...
                     1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1 ];

            if scalex/scaley>5 ||scaley/scalex>5
                warning('Scale factors differs significantly, mesh may of bad quality.');
            end
            obj.p  =  [scalex/3 0;0 scaley/2.5]*obj.p;
            while max(obj.triangleDiameters)>hmax
                obj.refineMesh;
            end
        end
    end

end

