function out = generate_random_edges_in_cylinder(N,container_radius,container_height,edge_length)
out = zeros(N,6);
for i = 1:N
    centroid = random_point_in_cylinder(container_radius,container_height);
    orientation = sphere_pick(1);
    first = centroid - edge_length*orientation/2;
    last = centroid + edge_length*orientation/2;
    out(i,:) = [first,last];
end

end