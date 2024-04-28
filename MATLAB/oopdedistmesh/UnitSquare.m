classdef UnitSquare < grid2D
    %UNITSQUARE UnitSquare class
    %   Creates an Unitsquare geometrie.

    % New Version fpr OCTAVEOOPDE: Simple hard codes initial grid.
    properties
        hmax
    end


    methods(Access = public)
        function obj = UnitSquare(hmax)
            % UnitSquare Constructor method
            %
            %     g = UnitSquare(hmax)
            %


            switch nargin
                case 0
                    % no empty call allowed
                    ME = MException('UNITSQUARE:MissingArgument',...
                               'Missing argument: Mesh-width h must be given.');
                    ME.throwAsCaller();
                case 1
                    if hmax <= 0
                        ME = MException('UNITSQUARE:MissformedArgument',...
                             'Missformed argument: Mesh-width h must be greater zero');
                        ME.throwAsCaller();
                    end
                    % Create mesh with 4 nodes and 2 elements by hand
                    obj.p = [0 0 1 1
                             0 1 0 1];

                    obj.e = [1     3     4     2
                             3     4     2     1
                             0     0     0     0
                             1     1     1     1
                             1     2     3     4];

                    obj.t = [3     2
                             2     3
                             1     4
                             1     1];

                    obj.hmax = hmax;

                    % Now refine it.
                    while obj.hmax > hmax
                      obj.refineMesh;
                    end
                otherwise
                    % empty
            end
        end
    end
end

