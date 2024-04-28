classdef Lshape < grid2D
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    
    methods
        function obj = Lshape(hmax)
            % Lshape
            % obj.lshape([hmax])
            % Creates an L shape domain with criss-cross
            % triangulation
            % (c) 2013 Uwe Prüfert
            
             switch nargin
                 case 0
                    % default constructor
                 case 1

                    obj.p = [0 0.5 1 0   0.5 1   0 0.5;...
                             0 0   0 0.5 0.5 0.5 1 1 ];

                    obj.e = [1   2   3 6 5 8 7   4;
                             2   3   6 5 8 7 4   1;...
                             0   0.5 0 0 0 0 0   0.5;...
                             0.5 1   1 1 1 1 0.5 1;...
                             1   1   2 3 4 5 6   6];

                    obj.t = [1 2 2 3 4 5 ;...
                             2 5 3 6 5 8;...
                             4 4 5 5 7 7;...
                             1 1 1 1 1 1];               

                    if ischar(hmax)
                        hmax = str2double(hmax);
                    end
                     
                    while max(obj.triangleDiameters) > hmax
                        obj.refineMesh;                        
                    end
                otherwise
                    obj.wrongNumberInputs.throwAsCaller;
            end
        end
    end    
end

