function s = compute_orientational_order2(edges)
% assert(size(edges,2) == 3)
% edges = edges./rwnorm(edges);

if size(edges,2) == 6
    temp = edges(:,1:3) - edges(:,4:6);
    temp = temp./rwnorm(temp);
    edges = temp;    
end

N = size(edges,1);

% orientational order
if N == 0
    aa = NaN(1,3);
    s = NaN;
elseif N > 0
    s = mean((3*edges(:,3).^2-1)/2);
end

end