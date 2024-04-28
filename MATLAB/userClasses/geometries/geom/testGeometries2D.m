% Test file for Circle class

% Define the circle using two intervals
x_range = [0 2*pi];
y_range = [-1 1];

% Create the Circle object
circle = Circle(1, x_range, y_range);

% Plot the circle
figure;
circle.plot;
axis equal;
title('Circle');

% Test the unit circle
circle = Circle();
figure;
circle.plot;
axis equal;
title('Unit Circle');

% Test a deformed circle
circle = Circle(1, x_range, y_range, 0.1);
figure;
circle.plot;
axis equal;
title('Deformed Circle');
