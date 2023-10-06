function two_point_edges = get_two_point_edges(labeled_edges)

insider_labels = unique(labeled_edges(:,7));

num_labels = numel(insider_labels);
two_point_edges = zeros(num_labels,6);

for i_lb = 1:num_labels
    e1 = labeled_edges(labeled_edges(:,7)==insider_labels(i_lb),1:3);
    e2 = labeled_edges(labeled_edges(:,7)==insider_labels(i_lb),4:6);
    
    p1 = e1(1,:);
    p2 = e2(end,:);
    
    two_point_edges(i_lb,1:6) = [p1,p2];
%     two_point_edges(i_lb,7) = insider_labels(i_lb);
end



end