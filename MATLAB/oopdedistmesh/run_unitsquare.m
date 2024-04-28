% Instantiate the UnitSquare object with hmax = 0.1
hmax = 0.1;
unitSquareObj = UnitSquare(hmax);

% Check if hmax property is set correctly
assert(unitSquareObj.hmax == hmax, 'hmax property is not set correctly');

% Validate the mesh properties
% Example: Check number of vertices (should be 4 for a unit square)
assert(size(unitSquareObj.p, 2) == 4, 'Number of vertices is incorrect');

% Additional validation checks can be performed on elements (t) and edges (e)
% based on the expected properties of a unit square mesh.

% If no assertion errors occur, the object is created correctly.
disp('UnitSquare object created successfully.');

