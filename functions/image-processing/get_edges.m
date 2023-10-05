function out = get_edges(rr_list)
N = numel(rr_list);
out = zeros(N,6);
for i = 1:N
    rr = rr_list{i};
    
    out(i,:) = [rr(1,:),rr(end,:)];
end

end