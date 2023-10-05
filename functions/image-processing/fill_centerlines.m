function new_cl_list = fill_centerlines(zstack,cl_list,d)

num_rods = numel(cl_list);
new_cl_list = cell(num_rods,1);
for i = 1:num_rods
    cl = cl_list{i};
    cl = rearrange_centerlines(cl,d);
    new_cl_list{i} = cl;

%     filled_rod = round(unique(cl,'rows','stable'));
%     filled_rod = cut_to_matrix_size(size(zstack),filled_rod);
% 
%     ind = sub2ind2(size(zstack),filled_rod);
%     true_ind = zstack(ind);
% 
%     ;
% 
%     new_cl_list{i} = filled_rod(true_ind,:);
end

end