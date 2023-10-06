function [TF,d_ij,i] = is_intersect2(edge,all_edge,rod_radius)
if isempty(all_edge)
    TF = false;
    return
end
TF = false;
N = size(all_edge,1);

for i = 1:N
    d_ij = distance_between_edges(edge,all_edge(i,:));
    if d_ij < 2*(1)*rod_radius
%         if d_ij < 2*(1+1e-6)*rod_radius
        TF = true;        
        return;
    end
end

end
