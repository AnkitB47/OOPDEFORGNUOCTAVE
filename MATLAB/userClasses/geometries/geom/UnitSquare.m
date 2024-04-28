classdef UnitSquare < grid2D
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
     
    
    methods(Access = public)
        function obj = UnitSquare(hmax)
            % UnitSquare
            % a version of unitsquare that creates criss-cross
            % Call
            %     g = UnitSquare(hmax)
            % Note that
            %
            %     g = UnitSquare()
            % creates an object with h = Inf.
            

            switch nargin
                case 0
                    % default constructor
                    obj.rectangle(0,1,0,1,inf);
                case 1
                    obj.rectangle(0,1,0,1,hmax);
                otherwise
                    %
            end
        end
    end    
end

