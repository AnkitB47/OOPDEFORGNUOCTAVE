% Define the distToSegment function
function d = distToSegment(p, v1, v2)
    % Calculate the distance from point p to the line segment defined by v1 and v2

    % Calculate the vector from v1 to v2 and from v1 to p
    v = v2 - v1;
    w = p - v1;

    % Calculate the projection of w onto v
    c1 = dot(w, v);
    c2 = dot(v, v);

    % If the projection falls outside the line segment, calculate the distance to the nearest endpoint
    if c1 <= 0
        d = norm(p - v1);
    elseif c2 <= c1
        d = norm(p - v2);
    else
        % Calculate the distance to the line segment
        b = c1 / c2;
        pb = v1 + b * v;
        d = norm(p - pb);
    end
end
