function inspect_rod_collisions(all_edges,rod_radius)

N = size(all_edges,1);
distance_matrix = Inf(N,N);
for i = 1:N
    edge_i = all_edges(i,:);
    
    for j = i+1:N
        edge_j = all_edges(j,:);
        if is_intersect2(edge_i,edge_j,rod_radius);
            fprintf('Intersection happens between %d and %d!\n',i,j)
            fprintf('Distance = %.2f\n',distance_between_edges(edge_i,edge_j))
            return
        end
    end
end

end