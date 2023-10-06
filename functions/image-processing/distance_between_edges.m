function [d,dist_vec,contact_site] = distance_between_edges(e1, e2)
% Each edge e1 and e2 is a 1x6 vector
% where the first three elements are the coordinates of the first point
% and the last three elements are the coordinates of the second point

% Define the edge points
p1 = e1(1:3);
q1 = e1(4:6);
p2 = e2(1:3);
q2 = e2(4:6);

% Direction vectors for the edges
d1 = q1 - p1;
d2 = q2 - p2;

% Vector from an endpoint of edge1 to an endpoint of edge2
r = p1 - p2;

a = dot(d1, d1);
b = dot(d1, d2);
c = dot(d2, d2);
e = dot(d1, r);
f = dot(d2, r);

s = 0.0; t = 0.0; % Parameters for parametric equations of two lines
det = a*c - b*b;  % Denominator in equations for s and t
s_numer = 0.0;    % Numerators in equations for s and t
t_numer = 0.0;

% Check if lines are not parallel
if abs(1.0 - b^2/a*c) > 1e-6
    s_numer = (b*f - c*e);
    t_numer = (a*f - b*e);
    
    % Check if s parameter for closest point on edge1 is within edge1
    if s_numer >= 0 && s_numer <= det
        s = s_numer / det;
    else
        if s_numer < 0
            s = 0;
        else
            s = 1;
        end
    end
    
    % Check if t parameter for closest point on edge2 is within edge2
    if t_numer >= 0 && t_numer <= det
        t = t_numer / det;
    else
        if t_numer < 0
            t = 0;
        else
            t = 1;
        end
    end
else
    s = -e/a;  % Comes from taking dot of e1 with a normal
    s = clip(s, 0.0, 1.0);
    
    t = (f + s * b) / c; % # Same as before
    t = clip(t, 0.0, 1.0);
end

% Distance vector
dist_vec = r + s*d1 - t*d2;

% Return the norm of the distance vector
d = norm(dist_vec);

contact_site = ((p1 + s*d1) + (p2 + t*d2))/2;

end