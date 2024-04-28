%% A boundary value problem in $R^1$
%
% Define a function. Since we want to use local functions, we cannot use a MATLAB script.
% Is it necessary to write a driver (main) function.
% Most of the code is identically with run_elliptic2D, exept market with (*).

function run_elliptic1D()

%% 
% Call the Elliptic constructor
    elliptic = Elliptic();
    
%%  
% (*) Assign an Interval object to grid property 
% Here, $I = (-1.1)$. 
% We use h = 0.05.
 

    elliptic.grid = Interval([-1,1],.05); 
   
%%
% (*) Assign Lagrange11D fem class with fem property  
    elliptic.fem = Lagrange11D();
    
%%
% Same as in 2D. Set boundary condition    
    elliptic.setBoundaryConditions('Dirichlet','0');

%% 
% (*)   Call initialize. To make the problem more interesting, we
% chose now parameters given by functions:  c = c(x), a = a(x) and source
% f = f(x). Use function-handles to give the parameters to initialize. 
    elliptic.initialize(@c,@a,@f);
    
%%
% Solve linear proeblem. Use 'LINEARGAUSS'.
% Try also 'LINEAR', and Algebraic
% multigrid solver by 'AMG'. 
% Note that AMG needs Ilupack installed.
    elliptic.solve('LINEARGAUSS');      
   
%%
% Plot the result. The option "LineStyle" =
% 'none' to supress printing the black lines is here senseless.
    elliptic.plot;
end
%% 
% (*) Coefficient function f and c defined by local functions. The source f is
% $f = 3$ if $x_1>0$, $x_2 >0$ and $x_3 >0$, otherwise $f=0$ in the "positive".
% The functions $c(x)$  and $a(x)$ impelements a Ball made of two components with
% different material properties.

function val = f(x1)
    val = 1*ones(1,length(x1));
    val(x1>0) = 0;   
end

function val = c(x1)
    val = 0.3*ones(1,length(x1));
    val(x1>0) = 3;   
end

function val = a(x1)
    val = 0.1*ones(1,length(x1));
    val(x1>0) = 1;
   
end