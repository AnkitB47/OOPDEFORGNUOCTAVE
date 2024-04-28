classdef Rectangle < grid2D
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods
        function obj = Rectangle(varargin)
        fprintf('GNU OCTAVE\n');
             switch nargin
                case 0
                    % default constructor
                otherwise
                    obj.rectangle(varargin{:});
            end
        end
    end
end

