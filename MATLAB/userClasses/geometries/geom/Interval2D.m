classdef Interval2D < handle
    properties
        x
        y
        hmax
    end

    methods
        function obj = Interval2D(x_range, y_range, hmax)
            if nargin < 3
                hmax = 0.1;
            end

            obj.x = linspace(x_range(1), x_range(2), ceil((x_range(2) - x_range(1)) / hmax) + 1);
            obj.y = linspace(y_range(1), y_range(2), ceil((y_range(2) - y_range(1)) / hmax) + 1);
            obj.hmax = hmax;
        end

        function plot(obj, f, varargin)
            if isempty(obj.x) || isempty(obj.y)
                error('Empty meshgrid X or Y');
            end

            % Evaluate the function f
            Z = f(obj.x, obj.y);

            % Check if Z is empty
            if isempty(Z)
                error('Empty function evaluation Z');
            end

            % Plot the surface
            figure;
            surf(obj.x, obj.y, Z, varargin{:});
            shading interp;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('Function Plot');
            colorbar;
        end
    end
end

