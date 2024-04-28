% Test file for Lagrange11D class

% Instantiate the Lagrange11D object with a specified mesh width
lagrange11D_object = Lagrange11D();

% Check if the object was created successfully
if isa(lagrange11D_object, 'Lagrange11D')
    disp('Lagrange11D object created successfully.');
else
    disp('Failed to create Lagrange11D object.');
end

