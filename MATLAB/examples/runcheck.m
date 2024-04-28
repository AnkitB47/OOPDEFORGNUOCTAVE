function runcheck()
    P = rand(100, 2); % generate a random set of 2D points
    dt = delaunay(P); % compute the Delaunay triangulation

    num_simplices = size(dt, 1); % get the number of simplices
    num_expected_simplices = round((size(P, 1) - 1) * size(P, 1) / 2); % expected number of simplices for a 2D Delaunay triangulation

    fprintf('Number of simplices: %d\n', num_simplices);
    fprintf('Expected number of simplices: %d\n', num_expected_simplices);

    if num_simplices ~= num_expected_simplices
        error("Unexpected number of simplices in Delaunay triangulation")
    end
end
