classdef Interval  < grid1D
     %Interval 1D FEM "Mesh"
     % Interval class implements an Interval I = [a,b] in means of p-e-t
     % structure.
     %
     % Arguments are
     % (i) (mandatory) a vector containing
     %    a) the start and end point of the interval or
     %    b) containing a discretization of the interval.
     % (ii) (optional) the maximal mesh size hmax.
     % Examples for calling (discretizing [a,b] = [0,1])
     %
     %    Interval([0, 1]);
     %    Interval([0, 0.2, 1]);
     %    Interval([0, 1], 0.1);
     %    Interval(linspace(0, 1, 100));
     %
     % In the presence of the second argument hmax, the mesh will be refined
     % by uniform refinement until the mesh size is smaller than hmax.
     %
     % This file is part of OOPDE for GNU Octave.




    methods(Access = public)
        function obj = Interval(varargin)
            switch nargin
                case 0
                   % Empty call. We will be here restrictive: Throw exception.
                   ME = MException("INTERVAL:EmptyCall",...
                        "Constructor must be called at least with one argument.");
                   ME.throwAsCaller;
                case 1
                    % Only boundary points or guess for discretization is given.
                    if ~isvector(varargin{1}) || length(varargin{1}) < 2
                        % Only one point, throw exception.
                        ME = MException("INTERVAL:MustBeVector",...
                             "First arguments must be a vector.");
                        ME.throwAsCaller;
                    else
                        % argument is a vector, we check if its well defined
                        varargin{1} = sort(varargin{1});
                        if max(varargin{1}) == min(varargin{1})
                            ME = MException("INTERVAL:MustBeMoreThanOnePoint",...
                                 "In [a,b], a must be different from b.");
                            ME.throwAsCaller;
                        end
                    end
                case 2
                    % [a,b] and hmax given.
                    % Test [a,b]. Right size?
                    if  ~isvector(varargin{1}) || length(varargin{1}) < 2
                        ME = MException("INTERVAL:MustBeVector",...
                             "First arguments must be a vector.");
                        ME.throwAsCaller;
                    else
                        % Argument is a vector, we check if it's well defined.
                        % Force it to be in the right order.
                        varargin{1} = sort(varargin{1});
                        if max(varargin{1}) == min(varargin{1})
                            ME = MException("INTERVAL:MustBeMoreThanOnePoint",...
                                 "In [a,b], a must be different from b.");
                            ME.throwAsCaller;
                        end
                    end
                    % Test hmax
                    if ~isscalar(varargin{2})
                        ME = MException("INTERVAL:MustBeScalar",...
                             "Second argument must be a positive scalar.");
                    elseif ~isnumeric(varargin{2})
                        ME = MException("INTERVAL:MustBeGreaterThanZero",...
                             "Argument hmax must be numeric.");
                        ME.throwAsCaller;
                    elseif isempty(varargin{2})
                        ME = MException("INTERVAL:MustBeGreaterThanZero",...
                             "Argument hmax must be non empty.");
                        ME.throwAsCaller;
                    elseif varargin{2} <= 0
                        ME = MException("INTERVAL:MustBeGreaterThanZero",...
                             "Argument hmax must be greater than zero.");
                        ME.throwAsCaller;
                    end
                otherwise
                    % We throw here an exception. Otherwise, obj.interval
                    % will do this later.
                    ME = MException("INTERVAL:TooManyArguments",...
                         "Too many arguments.");
                    ME.throwAsCaller;
            end
            % Call abstract interval method and give to it the entire
            %  parameter list.
            obj.interval(varargin{:});
        end
    end
end

