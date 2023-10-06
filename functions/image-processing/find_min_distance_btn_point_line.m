function [d_min,distance_vector,x_min] = find_min_distance_btn_point_line(point,rr)


dist = rwnorm(point - rr);
[quick_min,s_min] = min(dist);

s_end = size(rr,1);

vertices = rr(setdiff([s_min - 1, s_min, s_min + 1],[0,s_end+1]),:);
edges = diff(vertices);

d_min = Inf;
s_min1 = zeros(1,3);
s_min2 = zeros(1,3);
ee = zeros(1,3);

for ii = 1:size(edges,1)
    x0 = vertices(ii,:);
    x1 = vertices(ii+1,:);
    
    [candidate,dvec,xmin] = foo(point,x0,x1);
    if candidate < d_min
        d_min = candidate;
        distance_vector = dvec;
        x_min = xmin;        
    end
end

fig_on = 0;
if fig_on    
    close all;plot3v(rr,'o-');hold on;
    plot3v(point,'o','markersize',10);
end
% 
% if ~isequal(quick_min,d_min)
%     d_min
%     close all;plot3v(rr,'o-');hold on;
%     plot3v(point,'o','markersize',10);
%     ;
% end

end

function [d_min,distance_vector,xmin] = foo( point, x0, x1 )

s = dot( (point - x0), (x1 - x0) ) / norm(x1-x0)^2;

clip = @(x, low, high) max(low, min(x, high));
out_of_bounds = @(x, low, high) (x < low) | (x > high);

if s < 0
    xmin = x0;
elseif s > 1
    xmin = x1;
else
    xmin = x0 + s*(x1-x0);    
end
distance_vector = xmin - point;
d_min = norm(distance_vector);
end