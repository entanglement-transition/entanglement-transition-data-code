function [local_cell,nearby_rod_numberings]= find_nearby_rods(nosed_list,voxel_list,i_rod,thick_rad)
all_points = nosed_list(:,1:3);
noses = nosed_list(:,4);

rr_i = voxel_list{i_rod};
[cen_i,ori_i,slist_i] = get_line_coord(rr_i);

slist = sum( ( all_points - cen_i).*ori_i,2 );
dlist = rwnorm( all_points - ( cen_i + slist.*ori_i) );

length_estimate = max(slist_i) - min(slist_i);
indices_nearby_rods = dlist < thick_rad*5 & abs(slist) < length_estimate/2*1.5;

nearby_rod_numberings = unique(noses(indices_nearby_rods));
num_nearby_rods = numel(nearby_rod_numberings);
local_cell = voxel_list(nearby_rod_numberings);
% local_nosed_list = nosed_list(indices_nearby_rods,:);

% for j_rod = 1:num_nearby_rods
%     rod_number = nearby_rod_numberings(j_rod);
%     if rod_number == i_rod
%         continue;
%     end
%     
%     rr_j = voxel_list{rod_number};
%     
%     %     [~,d_ij] = find_min_distance_btn_discrete_lines(rr_i,rr_j);
% end

end