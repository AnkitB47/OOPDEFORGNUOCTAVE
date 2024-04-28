classdef ReactingFlow < pde & plotUtilsTimeDependent
    % class definiton for Buffoni Willcox example, a simple flame
    % simulation, cf. the paper "Projection-based model reduction for reacting flows"
    % by Marcelo Buffoni and Karen Willcox ...
    % It is a system of four diffusion-convection equations coupled by
    % a non-linear source term. The flow field is given as a static flow.
    
    properties(Dependent = true, GetAccess = public, SetAccess = private)
        temperature
        oxygen
        hydrogen
        water
    end
    
    
    
    % The "get methods"
    methods
        function val = get.temperature(obj)
            val = obj.y(1:obj.grid.nPoints,:);
        end
        function val = get.oxygen(obj)
            val = obj.y(1+obj.grid.nPoints:2*obj.grid.nPoints,:);
        end
        function val = get.hydrogen(obj)
            val = obj.y(1+2*obj.grid.nPoints:3*obj.grid.nPoints,:);
        end
        function val = get.water(obj)
            val = obj.y(1+3*obj.grid.nPoints:end,:);
        end
    end
    
    methods(Access = public)
        function initialize(obj,c,b)
            % Initialize the ractionFlow
            % We have only the diffusion coefficent and the flow
            % field, no source etc. 
           
            bcTempInlet = obj.grid.dirichletBC('min(950,max(300,sin(3*pi.*s-pi)*3500))');
            bcTempOverall = obj.grid.neumannBC('0');
            bcFuel = obj.grid.dirichletBC('max(0,sin(3*pi.*s-pi)*0.2)');
            bcOx = obj.grid.dirichletBC('max(0,sin(3*pi.*s-pi)*0.02)');
            bcProd = obj.grid.dirichletBC('1','0');
          

            bcNeumann = obj.grid.neumannBC('0'); 

            obj.grid.makeBoundaryMatrix(bcNeumann,bcNeumann,bcNeumann,bcFuel);
            [~,~,H,RF] = obj.fem.assemb(obj.grid);
           
            obj.grid.makeBoundaryMatrix(bcNeumann,bcNeumann,bcNeumann,bcOx);
            [~,~,~,RO] = obj.fem.assemb(obj.grid);
            
            obj.grid.makeBoundaryMatrix(bcNeumann,bcNeumann,bcNeumann,bcProd);
            [~,~,~,RP] = obj.fem.assemb(obj.grid);   
            
            obj.grid.makeBoundaryMatrix(bcTempOverall,bcTempOverall,...
                bcTempOverall,bcTempInlet);
            [~,~,~,RT] = obj.fem.assemb(obj.grid); 
             
            [K,~,~]= obj.fem.assema(obj.grid,c,0,0);
            C = obj.fem.convection(obj.grid,b);
            
            M = obj.mass();
            N = sparse(obj.grid.nPoints,obj.grid.nPoints);
            
            obj.A = -[K+C+1e8*(H'*H) N N N;
                      N K+C+1e8*(H'*H) N N;
                      N N K+C+1e8*(H'*H) N;
                      N N N K+C+1e8*(H'*H)];    
            
            obj.D = [ M N N N ;
                      N M N N;
                      N N M N;
                      N N N M];                  
            
            obj.b = [1e8*H'*RT;
                      1e8*H'*RF;
                      1e8*H'*RO;
                      1e8*H'*RP];                
                  
            obj.initialized = true;
        end
        
        
        
        function plot(obj,varargin)
            clf 
            subplot(2,2,1)     
            obj.grid.plot(obj.temperature(:,end),varargin{:});view(2)
            subplot(2,2,2)     
            obj.grid.plot(obj.oxygen(:,end),varargin{:});view(2)
            subplot(2,2,3)     
            obj.grid.plot(obj.hydrogen(:,end),varargin{:});view(2)
            subplot(2,2,4)     
            obj.grid.plot(obj.water(:,end),varargin{:});view(2)
        end   
        
        function movie(obj,varargin)
            try
                for k = 1:length(obj.time)                 
                    subplot(2,2,1)     
                    obj.grid.plot(obj.temperature(:,k),varargin{:});
                    view(0,90)
                    drawnow
                    subplot(2,2,3)     
                    obj.grid.plot(obj.oxygen(:,k),varargin{:});
                    view(0,90)
                    drawnow
                    subplot(2,2,2)     
                    obj.grid.plot(obj.hydrogen(:,k),varargin{:});
                    view(0,90)
                    drawnow
                    subplot(2,2,4)     
                    obj.grid.plot(obj.water(:,k),varargin{:});
                    view(0,90)
                    drawnow
                
                end
            catch ME
                fprintf(['A problem accured while playing the movie: ',...
                    ME.message,'\n']);
            end
        end
        
        function initialValue(obj)
            % solves the linear part of the system i.e. convection and
            % diffusion to obtain a good choice for an inital value
            % for the non linear solver.
            obj.y = -obj.A\obj.b;
        end
    end
    methods(Access = public,Static = true)
        % source term. For later use in DEIM it must be public...
        function d  = w(~,Y)
            % chemical source term for the Buffoni/Willcox paper project
            % w = g(Y,T)
            % Simple one step reaction is 2H2 + O2 => H2O
            % Fill the Data by using valus from cantera
            % possible is also a temperature dependent 
            % implementation for A
            % Molecularar weights: take this from Literature
            % Reaction data
            A = 5.5e11;
            E = 4.5e3;         
            nu = [2 1 -2]';
            wc= [2.016 31.9 18.0]';   
            R = 8.314472;
            rho = 1.39e-3;
            Q = 9800;
 
            
%             np = obj.grid.nPoints;
            % To use it 
            np = length(Y)/4;
            % this part is valid for all species, because of the use of a
            % global mechanism
            scalar_rate =...
                (rho*Y(np+1:2*np)/wc(1)).^nu(1).*...
                (rho*Y(2*np+1:3*np)/wc(2)).^nu(2).*A.*exp(-E./(R.*Y(1:np)));

            dY = -1/rho*[wc(1)*nu(1)*scalar_rate;...
                wc(2)*nu(2)*scalar_rate;...
                wc(3)*nu(3)*scalar_rate];   
            % for temperature is  dY_fuel/dt * Q
            dT = dY(2*np+1:3*np)*Q;   
            d = [dT;dY];
        end 
        
        
        function d  = dw(~,Y)
            % we need the Jacobian of our chemiscal source term, 
            % rather criptic...
            % the derivative of dw       
            A = 5.5e11;
            E = 4.5e3; 
            nu = [2 1 -2]';
            wc= [2.016 31.9 18.0]';   
            R = 8.314472;
            rho = 1.39e-3;
            Q = 9800;
            % np = obj.grid.nPoints;
            % Make it usualbe also for the scalar case
            np = length(Y)/4;
            % computed on paper and coded as a "matrix gradient"
            d_rate = [...
                    sparse(1:np,1:np, (rho/wc(1)*Y(np+1:2*np)).^(nu(1)).*...
                    (rho/wc(2)*Y(2*np+1:3*np)).^nu(2).*...
                    (A*exp(-E./(R.*Y(1:np)))).*(E./R./(Y(1:np).^2))),...
                    sparse(1:np,1:np, nu(1)*(rho/wc(1))*(rho*Y(np+1:2*np)/wc(1)).^(nu(1)-1).*...
                    (rho*Y(2*np+1:3*np)/wc(2)).^nu(2).*...
                    (A*exp(-E./(R.*Y(1:np))))),...                                
                    sparse(1:np,1:np, nu(2)*(rho/wc(2))*(rho/wc(1)*Y(np+1:2*np)).^(nu(1)).*...
                    (rho/wc(2)*Y(2*np+1:3*np)).^(nu(2)-1).*...
                    (A*exp(-E./(R.*Y(1:np))))),... 
                                         sparse(np,np)];
            % in 0D case this is a dyadic product...
            d =-1/rho* [Q*wc(3)*nu(3)*d_rate ;wc(1)*nu(1)*d_rate ;...
                    wc(2)*nu(2)*d_rate ;wc(3)*nu(3)*d_rate];
        end
    end
    
    
    
    
    methods(Access = protected)
        function val = df(obj,~,y)
            val = obj.A*y+obj.D*obj.w(0,y) + obj.b;   
        end
        
        
        function J = jacobian(obj, ~,y)          
            J = obj.A+obj.D*obj.dw(0,y); 
        end
        
    end
    
    
end