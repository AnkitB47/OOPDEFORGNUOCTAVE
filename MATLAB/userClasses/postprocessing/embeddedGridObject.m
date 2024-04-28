classdef embeddedGridObject < handle
    %% embededGridObject
    %
    % Class that implements submanifolds of grid objects.
    %
    % A submanifold is an object with dimension n-1 embedded in
    % $R^n$.

    %% Properties with SetAccess = private
    properties(SetAccess = private)
        data double
        manifold % The manifold (i.e. a grid object)
        % The submanifold is defined by their points, elements and their
        % dimension.
        p double % points of submanifold
        t double % t related with submanifold
        subdim double % dimension of submanifold
    end

    methods(Access = public)
        %% Public methods
        %
        function obj = embeddedGridObject(manifold,submanifold,varargin)
            %%
            % * embededGridObject IN:gridd,gridd[,double[,double]] OUT:self
            switch nargin
                case 0 % Empty object
                     % do nothing
                case 1
                    % Give through call or error..
                    if isa(manifold,'embededGridObject')
                        obj = manifold;
                    else
                        ME = MException('EMBEDEDGRIDOBJECT:WRONGNUMBERARGUMENST',...
                             'This call needs at least two arguments.');
                        ME.throwAsCaller();
                    end
                otherwise
                    % Two or more arguments
                    obj.manifold = manifold;
                    obj.subdim = submanifold.spaceDimension;
                    obj.t =  submanifold.t(1:end-1,:);
                    obj.p = submanifold.p;

                    for k = 1:obj.manifold.spaceDimension-submanifold.spaceDimension
                        if length(varargin)>=k
                            obj.p = [obj.p;varargin{k}];
                        else
                            obj.p = [obj.p;zeros(1,submanifold.nPoints)];
                        end
                    end
                    %
            end
        end

        function rotateSubmanifold(obj,varargin)
           %%
            % * rotateSubmanifold IN:self,cdouble[,double,double,double] OUT:self
            %
            % Method that rotates the points of a mesh. The method makes
            % only sense if the dimension of the object is 2 or 3. In 2D,
            % the argument is the rotation angle (positive means
            % counter-clockwise). In 3D, the first three arguments are
            % coordinates of the rotation axis (as vector from the origing)
            % and the fourth argument is the rotation angle.
            %
            % Call
            % 2D
            %
            %      rotateSubmanifold(alpha)
            %      rotateSubmanifold(n1,n2,n3,alpha)
            %
            % 3D
            %      rotateSubmanifold(x0,y0,z0,n1,n2,n3,alpha)
            %
            % Example
            %
            %       rotateSubmanifold(0,0,1,pi/4)
            %
            % rotates a 3D object around the z axis with angle
            % $\alpha = \pi/4$

            switch obj.manifold.spaceDimension
                case 2
                    switch length(varargin)
                        case 1 % only rotate around (0,0)
                            a = mod(varargin{3},2*pi);
                            RX = [cos(a) -sin(varargin{3})
                                  sin(a)  cos(varargin{3})];

                            obj.p = RX*obj.p;
                        case 3 % rotate around point (x0,y0)
                            % move object in (0,0)
                            obj.p = obj.p-[varargin{1}
                                           varargin{2}];


                            a = mod(varargin{3},2*pi);
                            % rotation matrix
                            RX = [cos(a) -sin(varargin{3})
                                  sin(a)  cos(varargin{3})];

                            % rotate object
                            obj.p = RX*obj.p+[varargin{1}
                                              varargin{2}];

                            % move back to (x0,y0)
                            obj.p = obj.p+[varargin{1}
                                           varargin{2}];
                        otherwise
                            obj.wrongNumberInputs.throwAsCaller();
                    end
                case 3
                    switch length(varargin)
                        case 4
                            % normalize vector
                            n1 = varargin{1}/norm([varargin{1},varargin{2},varargin{3}]);
                            n2 = varargin{2}/norm([varargin{1},varargin{2},varargin{3}]);
                            n3 = varargin{3}/norm([varargin{1},varargin{2},varargin{3}]);

                            a = mod(varargin{4},2*pi);
                            % rotation matrix
                            RN = [n1^2*(1-cos(a))+cos(a)        n1*n2*(1-cos(a))-n3*sin(a)  n1*n3*(1-cos(a))+n2*sin(a)
                                  n2*n1*(1-cos(a))+n3*sin(a)    n2^2*(1-cos(a))+cos(a)      n2*n3*(1-cos(a))-n1*sin(a)
                                  n3*n1*(1-cos(a))-n2*sin(a)    n3*n2*(1-cos(a))+n1*sin(a)  n3^2*(1-cos(a))+cos(a)];
                            % rotate
                            obj.p = RN*obj.p;
                        case 7
                            % bring it to (0,0,0)
                            obj.p = ob.p -[varargin{1}
                                           varargin{2}
                                           varargin{3}];

                            % normalize vector

                            n1 = varargin{4}/norm([varargin{4},varargin{5},varargin{6}]);
                            n2 = varargin{5}/norm([varargin{4},varargin{5},varargin{6}]);
                            n3 = varargin{6}/norm([varargin{4},varargin{5},varargin{6}]);

                            a = mod(varargin{7},2*pi);
                            % rotation matrix
                            RN = [n1^2*(1-cos(a))+cos(a)        n1*n2*(1-cos(a))-n3*sin(a)  n1*n3*(1-cos(a))+n2*sin(a)
                                  n2*n1*(1-cos(a))+n3*sin(a)    n2^2*(1-cos(a))+cos(a)      n2*n3*(1-cos(a))-n1*sin(a)
                                  n3*n1*(1-cos(a))-n2*sin(a)    n3*n2*(1-cos(a))+n1*sin(a)  n3^2*(1-cos(a))+cos(a)];
                            % rotate
                            obj.p = RN*obj.p+[varargin{1}
                                              varargin{2}
                                              varargin{3}];
                        otherwise

                    end
            end
        end

        function moveSubmanifold(obj,varargin)
            %%
            %    moveSubmanifold IN:self,double,[double,[double]];
            %
            % Moves submanifold. Arguments are x,y,z distance to move.


            switch obj.manifold.spaceDimension
                case 1
                    if isscalar(varargin{1})
                        obj.p(1,:) = obj.p(1,:)+varargin{1};
                    end
                case 2
                    if isscalar(varargin{1})&&...
                            isscalar(varargin{1})&&...
                            isscalar(varargin{1})
                        obj.p(1,:) = obj.p(1,:)+varargin{1};
                        obj.p(2,:) = obj.p(2,:)+varargin{2};

                    end
                case 3
                    if isscalar(varargin{1})&&...
                            isscalar(varargin{1})&&...
                            isscalar(varargin{1})
                        obj.p(1,:) = obj.p(1,:)+varargin{1};
                        obj.p(2,:) = obj.p(2,:)+varargin{2};
                        obj.p(3,:) = obj.p(3,:)+varargin{3};
                    end
            end
        end

        function setData(obj,data)
        %%
        % * setData IN:self,double OUT:none
        %
        % Uses data defined on manifold and interpolates it on submanifold.
        % If a point of the submanifold is not in the domain, its related
        % data is set to NaN.
        %


        obj.data = obj.manifold.interpolate(obj.manifold.p,obj.p,data);

        end

        function plot(obj,varargin)
            switch obj.manifold.spaceDimension
                case {2}
                    obj.manifold.plot();
                case 3
                    % 1D,2D 3/ eingebettet in 3D
                    obj.manifold.plotFaces([],'FaceColor','none');
                    switch obj.subdim
                        case 1
                            x = obj.p(1,:);
                            y = obj.p(2,:);
                            z = obj.p(3,:);

                            line('XData',x,...
                                'YData',y,...
                                'ZData',z,...
                                varargin{:})

                        case 2 % The sub is a 2D object embedded into IR^3.
                            x = obj.p(1,:);
                            y = obj.p(2,:);
                            z = obj.p(3,:);
                            if isempty(obj.data)
                                patch(...
                                    'faces',obj.t(1:3,:)',...
                                    'vertices',[x(:),y(:),z(:)],...
                                    'facecolor',[.51 .51 1],...
                                    'edgecolor',[0,0,0],...
                                    varargin{:});
                            else
                                patch(...
                                    'faces',obj.t(1:3,:)',...
                                    'vertices',[x(:),y(:),z(:)],...
                                    'facevertexcdata',obj.data(:),...
                                    'facecolor','interp',...
                                    'edgecolor',[0,0,0],...
                                    varargin{:});
                            end
                        case 3

                        otherwise
                    end
            otherwise
            end
        end

%         function plotSubmanifold(obj,varargin)
%             %%
%             % * plotSubmanifold IN:self[,char,char,...]
%             %
%             % Plots the submanifold. Arguments can be all arguments known
%             % by Matlab's patch class constructor method. Usuful can be
%             % used
%             %
%             % * 'FaceColor'
%             % * 'EdgeColor'
%             %
%             % Note that the submanifold will be displayed as it is, i.e. at
%             % the original position.
%             %
%             % This method can be used to visualize cuts through 3D objects
%             % easily. It uses some methods from gridd classes.
%             %
%
%
%             switch obj.subdim
%                 case 1
%                     %
%
%                     x = obj.p(1,:);
%                     y = obj.p(2,:);
%                     z = obj.p(3,:);
%
%                     if isempty(obj.data)
%                          line('XData',x,...
%                                 'YData',y,...
%                                 'ZData',z,...
%                                 varargin{:})
%                     else
%                          line('XData',x,...
%                                 'YData',y,...
%                                 'ZData',z,...
%                                 varargin{:})
%                     end
%
%                 case 2
%                     x = obj.p(1,:);
%                     y = obj.p(2,:);
%                     z = obj.p(3,:);
%                     if isempty(obj.data)
%                         % plots only the submanifold
%                         patch(...
%                             'faces',obj.t(1:3,:)',...
%                             'vertices',[x(:),y(:),z(:)],...
%                             'facecolor',[.51 .51 1],...
%                             'edgecolor',[0,0,0],...
%                             varargin{:});
%                     else
%                         % plots the submanifold and colored it with the
%                         % values of the data
%                         patch(...
%                             'faces',obj.t(1:3,:)',...
%                             'vertices',[x(:),y(:),z(:)],...
%                             'facevertexcdata',obj.data(:),...
%                             'facecolor','interp',...
%                             'edgecolor',[0,0,0],...
%                             varargin{:});
%                     end
%                 case 3
%                      % Not yet implemented
%                 otherwise
%             end
%         end
    end
end
%% Known issues
    %
    % Since setData method uses gridd.isPointInDomain method to set data to NaN
    % if a point from the submanifold is not in the domain, it can
    % be very slow.

