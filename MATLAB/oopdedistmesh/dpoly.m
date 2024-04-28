function d = dpoly(p, pv)
    % Copyright (C) 2024 Ankit Bhardwaj
    nvs = size(pv, 1) - 1;
    ds = zeros(size(p, 1), nvs);

    for iv = 1:nvs
        ds(:, iv) = donesegment(p, pv(iv:iv+1, :));
    end

    d = min(ds, [], 2);

    in_poly = inpolygon(p(:, 1), p(:, 2), pv(:, 1), pv(:, 2));
    d = (-1) .^ in_poly .* d;
end

function ds = donesegment(p, pv)
    v = diff(pv, 1);
    w = p - pv(1, :);

    c1 = sum(w .* v(1, :), 2);
    c2 = sum(v(1, :).^2);

    ds = sqrt(sum(w.^2, 2)) - abs(c1) ./ sqrt(c2);

    ds(c1 > 0 & c2 > c1) = sqrt(sum((p(c1 > 0 & c2 > c1, :) - pv(2, :)).^2, 2));
end

