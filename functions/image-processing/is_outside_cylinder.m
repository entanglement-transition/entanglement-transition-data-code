function TF = is_outside_cylinder(edge,container_radius,container_height)

tf1 = edge(3) < container_height;
tf2 = edge(3) > 0;
tf3 = edge(1)^2 + edge(2)^2 < container_radius^2;

tf4 = edge(6) < container_height;
tf5 = edge(6) > 0;
tf6 = edge(4)^2 + edge(5)^2 < container_radius^2;

TF =  tf1 & tf2 & tf3 & tf4 & tf5 & tf6;


end