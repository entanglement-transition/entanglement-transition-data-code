function [distance_vector,d_min,s_min1,s_min2,e1,e2] = find_min_distance_btn_discrete_lines( rr1,rr2 )

dmat = pdist2(rr1,rr2);
quick_min = min(dmat(dmat>0));
min_idx = find(dmat == quick_min);
[s_i,s_j] = ind2sub(size(dmat),min_idx);

max_i = size(rr1,1);
max_j = size(rr2,1);

vertices_i = rr1(setdiff([s_i - 1, s_i, s_i + 1],[0,max_i+1]),:);
vertices_j = rr2(setdiff([s_j - 1, s_j, s_j + 1],[0,max_j+1]),:);

edges_i = diff(vertices_i);
edges_j = diff(vertices_j);

d_min = Inf;
s_min1 = zeros(1,3);
s_min2 = zeros(1,3);
e1 = zeros(1,3);
e2 = zeros(1,3);
for i = 1:size(edges_i,1)
    e_i = [vertices_i(i,:) vertices_i(i+1,:)];
    for j = 1:size(edges_j,1)
        e_j = [vertices_j(j,:) vertices_j(j+1,:)];

        [dist_vec,so1,so2] = find_min_distance(e_i,e_j);        
%         [d,dvec,contact_point] = distance_between_edges(e_i,e_j);
        
        if norm(dist_vec) < d_min
            distance_vector = dist_vec;
            d_min = norm(dist_vec);
            s_min1 = so1;
            s_min2 = so2;   
            e1 = edges_i(i,:); e1 = e1/norm(e1);
            e2 = edges_j(j,:); e2 = e2/norm(e2);
        end
    end
end


end

