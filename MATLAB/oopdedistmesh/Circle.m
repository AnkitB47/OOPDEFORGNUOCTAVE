classdef Circle < grid2D
    %% Circle
    % Creates the mesh for a circle geometry
    % obj.circle() Unitcircle with hmax = 0.2
    % obj.circle(R) Circle with Radius = R and hmax = 0.2
    % obj.circle(R,hmax)
    % obj.circle(R,xshift,yshift) Circle with Radius = R,
    %    center (xshift, yshift) and hmax = 0.2
    % obj.circle(R,xshift,yshift,hmax)
    
    
    
    methods
        function obj = Circle(R,varargin)           
            switch nargin
                case 0
                     % default constructor
                case 1
                    obj.ellipse(R,R,R/10);
                case 3
                    obj.ellipse(R,R,varargin{:},R/10);
                case {2 4}
                    % Warn if hmax too large
                    if R < 2*varargin{end}
                        MException('CIRCLE:WRONGHMAX',...
                            ['The ratio R/hmax  = ',...
                            num2str(R/varargin{end}),...
                            '  < 2. Circle is not able to ',...
                            ' create deformed circles.']).throwAsCaller
                    end
                    obj.ellipse(R,R,varargin{:})%                 
                otherwise
                    obj.wrongNumberInputs.throw;
            end
        end
    end
end

