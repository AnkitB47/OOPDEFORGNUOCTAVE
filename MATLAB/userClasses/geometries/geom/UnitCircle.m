classdef UnitCircle < grid2D
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here



    methods
        function obj = UnitCircle(hmax)
            % obj.unitcircle
            % Uses drid2D.circle
            %
            % Call
            %   obj.unitcircle(hmax)
            %
            % Note that
            %     obj.UnitCircle()
            %
            % creates an empty object
            % (c) 2013 Uwe Prüfert

            switch nargin
                case 0
                    % default constructor
                case 1
                    obj.ellipse(1,1,0,0,hmax);
                otherwise
                    obj.wrongNumberInputs.throwAsCaller();
            end
        end
    end
end

