function d = fd_rectangle(p, x1, x2, y1, y2)
    % Signed distance function for a rectangle using drectangle
    d = -min(min(min(-y1 + p(:, 2), y2 - p(:, 2)), -x1 + p(:, 1)), x2 - p(:, 1));
end

