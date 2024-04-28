%% Durchbruchskurve
%
% Klassendefinition fuer das Model ...
%
%
% Die Differentialgleichung ist ein  System vom PDE-ODE Typ.
%
% PDE:
% 
% $$ \frac{\partial}{\partial t} c - D_{ax} \Delta c + \frac{U}{\epsilon_{sch}} \nabla  c + P_1 (c -\bar X)  = 0$$
% 
% mit Randbedingugnen 
% 
% $$c(0) = 1,\,\,\nabla c(1) = 0$$
%
% und Anfangswert $c(0,x) = 0$
%
% ODE
% 
% $$ \frac{\partial}{\partial t} \bar X + P2 (\bar X -c)= 0$$
%
% mit Anfangswert $\bar X(0,x) = 0$

classdef BreakThrough < pde & plotUtilsTimeDependent
    % Klasse fuer die Simulation der Durchbruchskurve.
    %%
    % Die Groessen $c$ und  $\bar X$ sind vom Loesungsvektor des
    % Gesamtsystems anhaengig. 
    properties(Dependent = true)
        c       % Konzentration
        X_bar   % mittlere Beladung
    end
    %%
    % Die get-Methoden fuer die abhaengigen Properties:
    methods
        function val = get.c(obj)
            val = obj.y(1:obj.grid.nPoints,:);
        end
        function val = get.X_bar(obj)
            val = obj.y(1+obj.grid.nPoints:end,:);
        end
    end
    
    %%
    % Immer gute Idee: Physikalische Daten als Properties zu definieren.
    properties(Access = public)
        D_ax = 8.89e-6
        U = 5.0
        eps_sch = 0.6
        D_eff = 2e-6
        Rho_s = 1500
        Rp = 0.005
        m = 1
    end
    %%
    % Diese beiden Methoden MUESSEN protected sein. 
    methods(Access = protected)
        function val = df(obj,~,y)
            % df ist hier die Standartformulierung fuer ein lineares System
            % von PDEs
            val = obj.A*y + obj.b;
        end
        
        function val = jacobian(obj,~,~)
            % jacobian ist hier einfach nur die Matirixvor dem linearen
            % Teil. Nicht notwendigerweise zu implementieren, macht den
            % Code aber viel schneller.
            val = obj.A;
        end
    end
    %%
    % Die oeffentlichen Methoden: initialize und solve.
    methods(Access = public)
        function initialize(obj)
            % In initialize MUSS obj.A und obj.b definiert werden.
            % Zusaetzlich sollte hier alles initialisiert werden, was hier
            % schon initialisiert werden kann: fem aendert sich nicht,
            % Randbedingungen aendern sich nicht mehr etc.
            
            
            % Damits es spaeter nicht so kompliziert wird, die Werte für
            % die Koeffizienten aus den physialischen Parametern ausrechnen.
            P2 = 15*obj.D_eff/obj.Rp/obj.Rp;              
            P1 = (1-obj.eps_sch)/obj.eps_sch*obj.Rho_s*P2; 
           
            
            % FEM  festlegen
            obj.fem = Lagrange11D;
            
            % Randbedingugnen festlegen: links Dirichlet, rechts Neumann.
            obj.grid.makeBoundaryMatrix(...
                obj.grid.dirichletBC( '1.0'),...
                obj.grid.neumannBC('0'));
            % Die ODE bekommt KEINE Randbedingung.
            
            % Die Matrizen assemblieren. Reihenfolge der Parameter:
            % obj.grid, Diffusion, Masse, Quelle, hier f=0.
            
            [K,M1,~] = obj.fem.assema(obj.grid,obj.D_ax,P1,0);
            
            % Nochmal fuer die   zweite Massematrx.
            [~,M2,~] = obj.fem.assema(obj.grid,0,P2,0);
            
            % Konvektionsterm. Parameter wieder obj.grid, Konvektion.
            C = obj.fem.convection(obj.grid,obj.U/obj.eps_sch);
            
            % Randmatrizen: Einziger Parameter: obj.grid. In obj.grid
            % steckt die Information zu den Randbedingugnen fuer jedes
            % Randsegment, hier zwei.
            [~,G,H,R] = obj.fem.assemb(obj.grid);
            
            % Der stif-spring Koeffizient für die Appriximation der
            % Dirichlet Radbedingung.
            lambda =  obj.fem.stiffSpring(K+M1);
            
            % Zum Fuellen der Systemmatrix benoetigen wir einen
            % Platzhalter, hier eine "Ultimate Sparse" Matrix 
            N = sparse(obj.grid.nPoints, obj.grid.nPoints);
            
            
            % Das PDE-ODE System formulieren: Erste Zeile PDE, zweite Zeile
            % ODE
            
            obj.A =  -[K+obj.m*M1+C+lambda*(H'*H)        -M1
                            -obj.m*M2                     M2];
                    
            obj.D = [obj.mass N
                      N obj.mass];
            
            % Kleiner Trick: G ist eignetlich hier fehl am Platze, 
            % aber das es ein Nullvektor pasender Groesse ist, setzen 
            % wir es hier ein      
            obj.b =   [lambda*H'*R
                          G];
             
            % Den Anfangswert fuer BEIDE Loesungskomponenten festlegen.          
            obj.y = zeros(2*obj.grid.nPoints,1);          
                      
            % Das "Bestaetigungsproperty" setzen.      
            obj.initialized = true;
        end
        %
        % Eigentlich nicht noetig, wir sichern uns aber gegen Fehlbenutzung
        % des Parameters "Loeser" ab. Erlaubt sind nur Loeser, die geeignet
        % sind.
        function solve(obj,arg)
            switch arg
                case {'ODE15S' 'EULERI' 'ODE23S' 'BDF2'}
                    solve@pde(obj,arg); 
                otherwise
                    MException('PDE:WRONGSOLVER',...
                        'Non suitable solver for this problem').throwAsCaller;
            end
        end   
    end
    
    %% check methods
    % 
    % For protecting the class for abuse, we restrickt to set fem only to
    % Lagrange11D and grid to Interval objects.
    % Its now not possible to use this class for 2D problems.
    methods(Static = true)
        function b = checkfem(fem)
            b = isa(fem,'Lagrange11D');
        end
        function b = checkgrid(grid)
            b = isa(grid,'Interval');
        end
    end
end

