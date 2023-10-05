function out = get_normal(x1,x2)
% get unit normal vector perpendicular to x1 and x2
out = cross(x1,x2,2)./rwnorm(cross(x1,x2,2));
end