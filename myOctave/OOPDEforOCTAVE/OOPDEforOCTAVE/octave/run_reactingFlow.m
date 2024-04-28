%% Reactiong flow
%
% Solves the "reacting flow" reactingflow stated in the paper
% "Projection-based model reduction for reacting flows"
% by Marcelo Buffoni and Karen Willcox 

function run_reactingFlow 
    %%
    % Call constructor
    reactingflow = ReactingFlow();
    %%
    % Create domain, mehs will be structured.
    reactingflow.grid = Rectangle(0,1.8,0,0.9,1); 
    %%
    % Refine the mesh uniformly
    reactingflow.grid.refineUniformly(5); 
    %%
    % Define $P^2$ elements.
    reactingflow.fem = Lagrange12D();
    %%
    %
    % If you want to try Lagrange 2 elements 
    % call    
    %         reactingflow.fem = lagrange22D();
    %         reactingflow.grid.extendMesh();    

    %%
    % Initialize the problem, Diffustiy $c=4$, $\vec b = [50,0]^\top$
    reactingflow.initialize(4,[50,0]');
    
    %%
    %
    reactingflow.initialValue();     

    %%
    %
    reactingflow.solve('STATIONARY');    
    
    %%
    %
    reactingflow.plot('EdgeColor','none') 
end