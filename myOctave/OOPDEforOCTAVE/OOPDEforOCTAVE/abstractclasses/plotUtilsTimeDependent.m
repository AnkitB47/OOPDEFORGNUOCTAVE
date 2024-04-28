%% plotUtilsTimeDependent
% Abstract class that defines additional methods for plotting time
% dependent problems. The user can add methods (by multy inheritance) to 
% its own user clases.
% 
% Example: 
%
% The class burgers implements the time-dependent Burgers equation. Add 
% plotUtilsTimeDependent to burgers:
%
%       classdef burgers < pde & plotUtilsTimeDependent
%
% Now you can use for instance the animation method inherited from 
% plotUtilsTimeDependent class:
%
%       burgers.animation('Command','grid on')
%
% Copy style = handle

%% Inheritance
%
% plotUtilsTimeDependent < handle

% plotTimeSpace plots solution over time-space using surf
% animation plots solution snapshots sequentially.

classdef (Abstract) plotUtilsTimeDependent < handle
  
    %#ok<*MCNPN>
     
    %% Public methods
    %
    methods
        function plotTimeSpace(obj,varargin) 
            %%
            %
            % * plotTimeSpace IN:self[,char]
            %
            % Plots the solution of a scalar time dependent problem with
            % spatial domain in 
            % $R^1$
            % as a surface plot over time and space.
            %
            % Arguments are all valid arguments that Matlabs surf method
            % understand and 'Command' or 'Commands' followed by char
            % arrays containing valid Matlab commands.
            %
            % Example:
            %
            %       burgers.plotTimeSpace('LineStyle','none',...
            %                             'Command','axis equal');
            %
            % plots the solution of burgers equation over the burgers.time
            % and burgers.grid.x. The edge lines are not printed and the
            % axis will be set to 'equal'. Further sencefull commands are
            % 
            % * view
            % 
            % * grid
            % 
             
            
            options = {'LineStyle','none'};
            commands = [];
            for k = 1:length(varargin)
                switch varargin{k}
                    case {'Commands' 'Command'}
                        k2 = 1;
                        for kl = k+1:length(varargin)                                                          
                            commands{k2} = varargin{kl};
                            k2 = k2 + 1;
                        end
                        break
                    otherwise
                        options{k+2} = varargin{k};
                end
            end        
            
            if size(obj.grid.p,1)>1 
                MException('plotUtilsTimeDependent:ONLY1D',...
                    'This method is only defined for 1D domains.').throwAsCaller;
            end
            [p,indx] = sort(obj.grid.p);             
            if any(size(obj.y)~=[length(p),length(obj.time)])
                MException('PLOTUTILSTIMEDEPENDENT:DIMENSION',...
                    ['Time and/or space discretization ',...
                    ' don''t match data.']).throwAsCaller;
            end
            surf(p,obj.time,obj.y(indx,:)',...
                options{:});
            for k2 = 1:length(commands)
                eval(commands{k2}); 
            end
            drawnow;
            xlabel('x','FontSize',12);
            ylabel('t','FontSize',12);
            zlabel('y(t,x)','FontSize',12);
        end
        
        function animation(obj,varargin) 
            %%
            %
            % * animation IN:self[,char]
            %
            % Plays a movie of the solution of a scalar time dependent 
            % problem 
            %
            % Arguments are all valid arguments that Matlabs surf method
            % understand and 'Command' or 'Commands' followed by char
            % arrays containing valid Matlab commands.
            %
            % Example:
            %
            %       problem.animation('-k','Command','grid on')
            %
            % Plays a movie of the solution of burgers equation. 
            % The color of the plot line was set to black and the grid was
            % switched on.
            % 
            commands = {};
            options = {};
            % Check options
            for k = 1:length(varargin)
                 switch varargin{k}
                     case {'Commands' 'Command'}
                         k2 = 1;
                         for kl = k+1:length(varargin)                                                          
                             commands{k2} = varargin{kl};
                             k2 = k2 + 1;
                         end
                         break
                     otherwise
                         options{k} = varargin{k};
                 end
            end                     
            
            dim = size(obj.grid.p,1);
            switch dim
                case 1 % animation of line plots
                    xmin = min(obj.grid.p);
                    xmax = max(obj.grid.p);
                    ymax = max(max(obj.y));
                    ymin = min(min(obj.y));           
                    for k = 1:length(obj.time)  
                        cla
                        obj.grid.plot(obj.y(:,k),options{:});
                        axis([xmin xmax ymin ymax]);
                        xlabel('x','FontSize',12);
                        ylabel('y(t,x)','FontSize',12);
                        title(['t = ',num2str(obj.time(k),'%3.2f')],...
                            'FontWeight','normal',...
                            'FontSize',12);
                        for k2 = 1:length(commands)
                            eval(commands{k2}); 
                        end
                        drawnow;
                    end                    
                case 2 % animation of surface plots         
                    xmin = min(obj.grid.p(1,:));
                    xmax = max(obj.grid.p(1,:));
                    ymax = max(obj.grid.p(2,:));
                    ymin = min(obj.grid.p(2,:)); 
                    zmax = max(max(obj.y));
                    zmin = min(min(obj.y));                    
                    for k = 1:length(obj.time)  
                        try
                            
                            maxValue = max(max(obj.y));
                            minValue = min(min(obj.y))-5;
                            set(gca,'Clim',[minValue, maxValue]);
                            colorbar('off')
                            axis([xmin xmax ymin ymax zmin zmax]);
                            obj.grid.plot(obj.y(:,k),options{:});
                            xlabel('x_{1}','FontSize',12);
                            ylabel('x_{2}','FontSize',12);
                            zlabel('y(t,x)','FontSize',12);
                            title(['t = ',num2str(obj.time(k),'%3.2f')],...
                                'FontWeight','normal',...
                                'FontSize',12) 
                            for k2 = 1:length(commands)                             
                                eval(commands{k2});
                            end
                            drawnow;
                        catch ME
                            fprintf(['A problem while playing occurs: ',...
                                ME.message,'\n']);
                        end
                    end
                case 3 % Animation of an outer faces plot
                    xmin = min(obj.grid.p(1,:));
                    xmax = max(obj.grid.p(1,:));
                    ymax = max(obj.grid.p(2,:));
                    ymin = min(obj.grid.p(2,:)); 
                    zmax = max(obj.grid.p(3,:));
                    zmin = min(obj.grid.p(3,:));
                    axis([xmin xmax ymin ymax zmin zmax]);    
                    for k = 1:length(obj.time)    
                        xlabel('x_{1}','FontSize',12);
                        ylabel('x_{2}','FontSize',12);
                        zlabel('x_{3}','FontSize',12);
                        title(['t = ',num2str(obj.time(k),'%3.2f')],...
                            'FontWeight','normal',...
                            'FontSize',12)
                        obj.grid.plotFaces(obj.y(:,k),options{:}); drawnow;
                        maxValue = max(max(obj.y));
                        minValue = min(min(obj.y))-5;
                        set(gca,'Clim',[minValue, maxValue]);
                        for k2 = 1:length(commands)
                            eval(commands{k2})
                        end
                    end                    
                otherwise
                    MException('plotUtilsTimeDependent:WRONGDIMENSION',...
                        'Domain must have dimension 1, 2, or 3.').throwAsCaller
            end
        end
    end    
end

