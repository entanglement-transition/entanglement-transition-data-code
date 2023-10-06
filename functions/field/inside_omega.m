function [omega,num_labels] = inside_omega(centerlines,cx,cy,cz,R,h)
% R: radius of bounding spheres
% h: grid spacing

num_rods = length(centerlines);
temp_cell = cell(size(centerlines));
for i = 1:num_rods
    cl = centerlines{i};
    %     rr = unique(rr,'rows','stable');
    N_cl = size(cl,1);
    temp_cell{i} = [cl(1:end-1,:),cl(2:end,:),i*ones(N_cl-1,1)];
end
labeled_centerline_edges = vertcat(temp_cell{:});
clear temp_cell

center = [cx,cy,cz];

I_cl = rwnorm(labeled_centerline_edges(:,1:3) - center) < R;
labeled_edges_bound = labeled_centerline_edges(I_cl,:);

insider_labels = labeled_centerline_edges(I_cl,7);
insider_labels = unique(insider_labels);

num_labels = numel(insider_labels);
omega = labeled_edges_bound;

end