% Define the rectangle size and channel dimensions
rect_width = 4;
rect_height = 8;
channel_width = rect_width;
channel_height = rect_height / 2;

% Calculate the coordinates of the rectangle vertices
rect_x = [0, rect_width, rect_width, 0, 0];
rect_y = [0, 0, rect_height, rect_height, 0];

% Calculate the coordinates of the channel vertices
channel_x = [rect_width/2 - channel_width/2, rect_width/2 - channel_width/2, rect_width/2 + channel_width/2, rect_width/2 + channel_width/2, rect_width/2 - channel_width/2];
channel_y = [0, channel_height, channel_height, 0, 0];

% Add origin for the channel
channel_x = [channel_x, rect_width/2];
channel_y = [channel_y, 0];

% Create rectangle and channel arrays for patch
rect_array = [rect_x, rect_y];
channel_array = [channel_x, channel_y];

% Create rectangle and channel patches
rect_patch = patch(rect_array(:,1), rect_array(:,2), [0.8, 0.8, 0.8], "FaceAlpha", 0.5);
channel_patch = patch(channel_array(:,1), channel_array(:,2), [0.3, 0.7, 0.3], "FaceAlpha", 0.5);

% Translate the channel to be below the rectangle
move_rect_channel = [rect_patch, channel_patch];
move_rect_channel = movpobj(move_rect_channel, [0, -rect_height]);

% Plot the result
figure();
hold on;
drawnow;

% Remove patch outlines
set(rect_patch, "EdgeColor", "none");
set(channel_patch, "EdgeColor", "none");

% Plot the final patches
plot(move_rect_channel);

% Add a legend
legend("Rectangle", "Channel");

% Set axis limits
axis([-1, rect_width + 1, -1, rect_height + channel_height + 1]);

% Set aspect ratio
axis equal;

hold off();
