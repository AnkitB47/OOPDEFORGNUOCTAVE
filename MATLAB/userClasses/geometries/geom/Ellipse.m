%% ellipse
%
%%
% Class that meshes an ellipse
%
% ellipse(A,B,x0,y0,hmax)
% 
%%
% $$(x-x_0)^2/A+(y-y_0)^2/B = 1$$
%

classdef Ellipse < grid2D 
    
    methods(Access = public)
        function obj = Ellipse(varargin) 
            switch nargin
                case 0
                    % default constructor
                otherwise
                    obj.ellipse(varargin{:});
            end
        end
    end    
end

