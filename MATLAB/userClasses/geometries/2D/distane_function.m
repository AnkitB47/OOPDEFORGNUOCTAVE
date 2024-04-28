% Distance function for a quadrant circle
function fd = distance_function(p)
    % Quadrant circle with radius 1
    xc = 0;
    yc = 0;
    r = 1;

    % Calculate distance from each point to the quadrant circle
    d = sqrt((p(:, 1) - xc).^2 + (p(:, 2) - yc).^2) - r;
end
