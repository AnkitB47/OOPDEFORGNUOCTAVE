% Test file for Lagrange12D class

% Instantiate the Lagrange12D object with a specified mesh width
hmax = 0.05; % Define your desired mesh width here
lagrange12D_object = Lagrange12D(hmax);

% Check if the object was created successfully
if isa(lagrange12D_object, 'Lagrange12D')
    disp('Lagrange12D object created successfully.');
else
    disp('Failed to create Lagrange12D object.');
end

