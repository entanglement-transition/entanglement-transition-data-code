function [distance_vector,so1,so2] = find_min_distance( edge1, edge2 )

if isempty(edge2)
    distance_vector = NaN;
    return
end

x1 = edge1(:,1:3); e1 = edge1(:,4:6) - x1;
x2 = edge2(:,1:3); e2 = edge2(:,4:6) - x2;

x1 = repmat(x1,[size(x2,1) 1]);
e1 = repmat(e1,[size(x2,1) 1]);

e1e1 = dot(e1, e1,2);
e1e2 = dot(e1, e2,2);
e2e2 = dot(e2, e2,2);

x1e1 = dot(x1, e1,2);
x1e2 = dot(x1, e2,2);
x2e1 = dot(e1, x2,2);
x2e2 = dot(x2, e2,2);

s = 0.0;
t = 0.0;

clip = @(x, low, high) max(low, min(x, high));
out_of_bounds = @(x, low, high) (x < low) | (x > high);

parallel = abs(1.0 - e1e2.^2 / (e1e1 .* e2e2)) < 1e-8;
if parallel
    % Some are parallel, so do processing
    t = (x2e1 - x1e1) / e1e1;  % Comes from taking dot of e1 with a normal
    t = clip(t, 0.0, 1.0);
    s = (x1e2 + t * e1e2 - x2e2) / e2e2;  % Same as before
    s = clip(s, 0.0, 1.0);
else
    % Using the Cauchy-Binet formula on eq(7) in docstring referenc
    s = (e1e1 .* (x1e2 - x2e2) + e1e2 .* (x2e1 - x1e1)) ./ (e1e1 * e2e2 - (e1e2) .^ 2);
    t = (e1e2 .* s + x2e1 - x1e1) ./ e1e1;
    
    if out_of_bounds(s, 0.0, 1.0) | out_of_bounds(t, 0.0, 1.0)
        %             # potential_s = -100.0
        %             # potential_t = -100.0
        %             # potential_d = -100.0
        %             # overall_minimum_distance = 1e20
        %             # Fill in the possibilities
        potential_t = (x2e1 - x1e1) ./ e1e1;
        s = 0.0;
        t = clip(potential_t, 0.0, 1.0);
        potential_d = rwnorm(x1 + e1 * t - x2);
        overall_minimum_distance = potential_d;
        
        potential_t = (x2e1 + e1e2 - x1e1) / e1e1;
        potential_t = clip(potential_t, 0.0, 1.0);
        potential_d = rwnorm(x1 + e1 * potential_t - x2 - e2);
        if potential_d < overall_minimum_distance
            s = 1.0;
            t = potential_t;
            overall_minimum_distance = potential_d;
        end
        
        potential_s = (x1e2 - x2e2) / e2e2;
        potential_s = clip(potential_s, 0.0, 1.0);
        potential_d = rwnorm(x2 + potential_s .* e2 - x1);
        if potential_d < overall_minimum_distance
            s = potential_s;
            t = 0.0;
            overall_minimum_distance = potential_d;
        end
        
        potential_s = (x1e2 + e1e2 - x2e2) ./ e2e2;
        potential_s = clip(potential_s, 0.0, 1.0);
        potential_d = rwnorm(x2 + potential_s .* e2 - x1 - e1);
        if potential_d < overall_minimum_distance
            s = potential_s;
            t = 1.0;
        end
    end
end

so1 = x1 + t * e1;
so2 = x2 + s * e2;

distance_vector = x2 + s * e2 - x1 - t * e1;


end