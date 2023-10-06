function all_edges = generate_intersecting_random_rods_in_cylinder(num_rods,container_radius,container_height,edge_length)

all_edges = [];
iterations = 0;

while 1
    edge = generate_random_edges_in_cylinder(1,container_radius,container_height,edge_length);
    distance_lower_bound = extended_line_distances(edge,all_edges);
%     neighbor_edges = all_edges(distance_lower_bound < 10*rod_radius,:);

    
    if ~is_outside_cylinder(edge,container_radius,container_height)
        continue
    else
        all_edges(end+1,:) = edge;
    end
    
    iterations = iterations + 1;
    if iterations >= num_rods
        break;
    end
end
% inspect_rod_collisions(all_edges,rod_radius);

end