%% ellipticAdaptive
% Class definition for diffusion equation. 
% It proviedes a solver with an adaptive mesh refinement,
% solveAdaptive.
% For that, we derive it from linearDiffusion class and add three public
% properties to store the coefficient functions and add the an adaptive 
% solver method.
%% Copy style
% handle

%% Inheritance
% ellipticAdaptive < linearDiffusion  
classdef EllipticAdaptive < Elliptic
    
    properties(Access = public)
        
        %% Properties with Access = public        
        %
        % * c, a, f (double|char|function_handel)
        %
        % For the adaptive solver, we need properties to store the 
        % coefficient (functions).
        c 
        a 
        f 
    end  
    
    methods ( Access = public )
        function initialize(obj,varargin)
            %%
            % * initialize
            % IN:self,double|char|function_handel,
            % double|char|function_handel,
            % double|char|function_handel
            % 
            % Intitializes the object. Call it w/o or with three arguments.
            % When calling w/o arguments, the properties obj.c obj.a and obj.f must be set
            % before.
            % 
            % Example:
            %%
            % 
            % 
            % (i) Predefined coefficent functions
            %             
            %
            %       obj.c = '1';
            %       obj.a = 0;
            %       obj.f = @f;
            %       obj.initialize()
            % 
            % (ii) Set coefficient function via initialize method
            %
            %       obj.initialize('1',0,@f)
            
            % (C) 2014 by Uwe Prüfert
            switch nargin
                case 4
                    obj.c = varargin{1};
                    obj.a = varargin{2};
                    obj.f = varargin{3};
                case 1
                    if (isempty(obj.c)||isempty(obj.a)||isempty(obj.f))                       
                        MException('PDE:COEFFIENCETS',...
                            'No coefficient functions set.').throwAsCaller
                    end
                otherwise
                    MException('PDE:WRONGNUMBERARGUMENTS',...
                        'The number of arguments must be zero or three.').throwAsCaller
            end
            b = zeros (obj.grid.spaceDimension ,1) ;        
            initialize@pde(obj,1,obj.c,b,obj.a,obj.f) ;
            l = obj.fem.stiffSpring(obj.K + obj.M );            
            obj.A = -(obj.K + obj.M + l*(obj.H'* obj.H ) + obj.Q ) ;
            obj.b = obj.F + l *(obj.H'* obj.R) + obj.G ;
        end
    
        function solveAdaptive(obj,tol) 
            %%
            % * solveAdaptive IN:self,double OUT:self
            %
            % Solves the problem by using an adaptive mesh refinement.
            % 
            % The number of refinements can be cotrolled by setting
            % solverOptions.maxit, but the maximum number of iterations 
            % is 10.
            %
            % Example
            %
            %       heatTransfer.solveAdaptive(); 
            %
            %       heatTransfer.solveAdaptive(1e-2); 
            
            % (c) 2015 by Uwe Prüfert            
            
            maxit = min(10,obj.solverOptions.maxit);
            niter = 1;  
            if ~obj.initialized
                MException('PDE:NOTINITIALIZED',...
                    'The object is not initialized.').throwAsCaller
            end
            while true 
                 
                obj.y = zeros(obj.grid.nPoints,1);
                obj.solve() ;
                % make it dimension independent                
                switch class(obj.c)
                    case 'function_handle'                        
                        switch obj.grid.spaceDimension
                            case 1 % 1D domain
                                cc = obj.c(obj.grid.p(1,:));  
                            case 2
                                 cc = obj.c(obj.grid.p(1,:),...
                                     obj.grid.p(2,:));   
                            case 3
                                 cc = obj.c(obj.grid.p(1,:),...
                                     obj.grid.p(2,:),obj.grid.p(2,:));   
                            otherwise
                        end
                    case 'double' % dimension independent
                        if  isscalar(obj.c)
                            cc = obj.c*ones(obj.grid.nPoints,1);                            
                        else
                            MException('PDE:WRONGCLASS',...
                                ['For adaption, all data must be',...
                                 ' function_handles or scalar']).throw ;
                        end
                    case 'char'
                        try
                            x = obj.grid.x;
                            y = obj.grid.y;
                            z = obj.grid.z;
                        catch
                            % don't catch, try 'n' error code...
                        end
                        cc = eval(obj.c);
                        if  isscalar(obj.c)
                            cc = cc*ones(size(obj.grid.p,2),1); 
                        elseif length(cc)==obj.grid.nPoints
                            % okay
                        else
                            MException('PDE:WRONGCLASS',...
                                ['For adaption, all data must be',...
                                 ' function_handles or scalar']).throw ;
                        end
                        
                    otherwise
                        MException('PDE:WRONGCLASS',...
                            ['For adaption, all data must be',...
                             ' function_handles or double']).throw ;
                end 
                switch class(obj.a)
                    case 'function_handle'
                        switch obj.grid.spaceDimension
                            case 1 % 1D domain
                                aa = obj.a(obj.grid.p(1,:));  
                            case 2
                                aa = obj.a(obj.grid.p(1,:),...
                                     obj.grid.p(2,:));   
                            case 3
                                aa = obj.a(obj.grid.p(1,:),...
                                     obj.grid.p(2,:),obj.grid.p(2,:));   
                            otherwise
                        end
                    case 'double'
                        if  isscalar(obj.a)
                            aa = obj.a*ones(obj.grid.nPoints,1);
                        else
                            MException('PDE:WRONGCLASS',...
                            ['For adaption, all data must be',...
                             ' function_handles or scalar']).throw ;
                        end
                    case 'char'
                        try
                            x = obj.grid.x;
                            y = obj.grid.y;
                            z = obj.grid.z;
                        catch
                            % don't catch, try 'n' error code...
                        end
                        aa = eval(obj.a);
                        if  isscalar(aa)
                            aa = aa*ones(obj.grid.nPoints,1); 
                        elseif length(aa)==obj.grid.nPoints
                            % okay
                        else
                            MException('PDE:WRONGCLASS',...
                                ['For adaption, all data must be',...
                                 ' function_handles or scalar']).throw ;
                        end
                    otherwise
                        MException('PDE:WRONGCLASS',...
                            ['For adaption, all data must be',...
                            ' function_handles or double']).throw ;
                end 
                switch class(obj.f)
                    case 'function_handle'
                        switch obj.grid.spaceDimension
                            case 1 % 1D domain
                                ff = obj.f(obj.grid.p(1,:));  
                            case 2
                                ff = obj.f(obj.grid.p(1,:),...
                                     obj.grid.p(2,:));   
                            case 3
                                ff = obj.f(obj.grid.p(1,:),...
                                     obj.grid.p(2,:),obj.grid.p(3,:));   
                            otherwise
                        end
                    case 'double'
                        if  isscalar(obj.f)
                            ff = obj.a*ones(obj.grid.nPoints,1);
                        else
                            MException('PDE:WRONGCLASS',...
                            ['For adaption, all data must be',...
                             ' function_handles or scalar']).throw ;
                        end
                    otherwise
                        MException('PDE:WRONGCLASS',...
                            ['For adaption, all data must be',...                      
                             ' function_handles or doubel']).throw ;
                end  

                E =  obj.errorInd(cc,aa,ff,1,1,0);                
                if max(E)<tol || niter>maxit
                    break
                end
                indx = obj.selectElements2Refine(E,0.125);
                obj.grid.refineMesh(indx); 
                obj.initialize();
                niter = niter + 1;             
                fprintf(['Solve on adaptive refined mesh with ',...
                    num2str(obj.grid.nPoints),' nodes.\n']);
            end           
        end
    end
end
