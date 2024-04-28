% Size function for a non-uniform triangular mesh
function fh = size_function(p)
    % Compute size based on distance from the center of the circle
    r = sqrt(p(:, 1).^2 + p(:, 2).^2);

    % Adjust size function to promote non-uniform mesh
    h = 0.05 + 0.3 * abs(r - 1);
end
