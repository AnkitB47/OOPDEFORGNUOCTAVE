classdef CustomRectangle < grid2D
    properties
        width
        height
    end

    methods
        function obj = CustomRectangle(varargin)
            fprintf('GNU OCTAVE\n');
            switch nargin
                case 0
                    % default constructor
                otherwise
                    obj = obj@grid2D(varargin{:});
            end
        end
        % Public method to get the fd and fh functions
        function [fd, fh] = getFunctions(obj)
            fd = @obj.fd;
            fh = @obj.fh;
        end
    end

    methods(Access = protected)
        function d = fd(obj, x, y)
            % Signed distance function for the custom rectangle
            width = obj.width;
            height = obj.height;

            % Compute distances to the sides of the custom rectangle
            dx = max(abs(x) - width / 2, 0);
            dy = max(abs(y) - height / 2, 0);

            % Combine distances using the max function to create the union
            d = sqrt(dx.^2 + dy.^2);
        end

        function h = fh(obj, x, y)
            % Mesh density function (adjust as needed)
            h = 0.1;  % Desired element size
        end
    end

     methods(Access = public)  % Change access to public
        function vertices = getVertices(obj)
            % Get the vertices of the custom rectangle
            width = obj.width;
            height = obj.height;

            % Vertices of the custom rectangle
            vertices = [-width / 2, -height / 2;
                         width / 2, -height / 2;
                         width / 2, height / 2;
                        -width / 2, height / 2];
        end
    end
end
