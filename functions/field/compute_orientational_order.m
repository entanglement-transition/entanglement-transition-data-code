function s = compute_orientational_order(edges)
assert(size(edges,2) == 6)
temp = edges(:,1:3) - edges(:,4:6);
temp = temp./rwnorm(temp);
orientations = temp;    

N = size(orientations,1);

% orientational order
if N == 0
    aa = NaN(1,3);
    s = 0;
elseif N > 0
    % edges0 = edges(labels,:);
    Q = zeros(3,3);
    for k = 1:size(orientations,1)
        v = orientations(k,:);
        Q = Q + kron(v,v');
    end
    Q2 = (3*Q/N - eye(3))/2;
    [V,D] = eig(Q2);
    
    [~,I] = max((diag(D)));    
%     [~,I] = max(abs((diag(D))));
    
    s = D(I,I);
end

end