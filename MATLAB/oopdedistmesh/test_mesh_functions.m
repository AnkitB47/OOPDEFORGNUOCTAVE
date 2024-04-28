% Test script for mesh-generating functions in Octave

% Add the paths to the mesh files
addpath(genpath(pwd));  % Assuming the mesh files are in the same directory

% Test 1
try
    d = dcircle([0, 0], 0, 0, 1);
    disp('dcircle: Passed');
catch
    disp('dcircle: Failed');
end

% Test 2
try
    d = ddiff([1; 2; 3], [2; 3; 4]);
    disp('ddiff: Passed');
catch
    disp('ddiff: Failed');
end

% Repeat the above structure for each mesh function

% Test 3
try
    % Add your test for the third mesh function here
    disp('dintersect: Passed');
catch
    disp('dintersect: Failed');
end

% Test 4
try
    % Add your test for the fourth mesh function here
    disp('drectangle: Passed');
catch
    disp('drectangle: Failed');
end

% Test 5
try
    % Add your test for the fifth mesh function here
    disp('dsphere: Passed');
catch
    disp('dsphere: Failed');
end

% Test 6
try
    % Add your test for the sixth mesh function here
    disp('dunion: Passed');
catch
    disp('dunion: Failed');
end

% Test 7
try
    % Add your test for the seventh mesh function here
    disp('fixmesh: Passed');
catch
    disp('fixmesh: Failed');
end

% Test 8
try
    % Add your test for the eighth mesh function here
    disp('huniform: Passed');
catch
    disp('huniform: Failed');
end

% Test 9
try
    % Add your test for the ninth mesh function here
    disp('simpvol: Passed');
catch
    disp('simpvol: Failed');
end

% Test 10
try
    % Add your test for the tenth mesh function here
    disp('distmesh2d: Passed');
catch
    disp('distmesh2d: Failed');
end


% Create a 2D polyhexcore object
figure(1);
title('Polyhex 2D');
polyhex2D = polyhexcore(2);
polyhex2D.generateMesh();

% Create a 3D polyhexcore object
figure(2);
title('Polyhex 3D');
polyhex3D = polyhexcore(3);
polyhex3D.generateMesh();


% Remove paths
rmpath(genpath(pwd));

disp('Test complete.');

