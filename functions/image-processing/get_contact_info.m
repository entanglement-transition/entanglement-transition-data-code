function contact_table = get_contact_info(voxel_list,stack,radius_threshold)

num_tracked = numel(voxel_list);

distance_lowerbound = Inf(num_tracked,num_tracked);

nosed_list = segment_cells_to_nosed_list(voxel_list);
all_points = nosed_list(:,1:3);
noses = nosed_list(:,4);


i = [];
j = [];
p_i1 = [];
p_i2 = [];
s_i = [];
s_j = [];
Q_i = [];
Q_j = [];
%%
num_contact_voxels = [];

t_start = tic;
for i_rod = 1:num_tracked
    rr_i = voxel_list{i_rod};
    [cen_i,ori_i,slist_i] = get_line_coord(rr_i);
    fake_line_i = [cen_i,cen_i+10*ori_i];
    
    slist = sum( ( all_points - cen_i).*ori_i,2 );
    dlist = rwnorm( all_points - ( cen_i + slist.*ori_i) );
    length_estimate = max(slist_i) - min(slist_i);
    indices_nearby_rods = dlist < radius_threshold*5 & abs(slist) < length_estimate/2*1.5;
    nearby_rod_numberings = unique(noses(indices_nearby_rods));
    num_nearby_rods = numel(nearby_rod_numberings);
    local_nosed_list = nosed_list(indices_nearby_rods,:);
    
    [rod_i,ind_i] = get_rod_from_centerline(rr_i,radius_threshold,stack);
    
    contacts_i = [];
    for j_rod = 1:num_nearby_rods
        rod_number = nearby_rod_numberings(j_rod);
        if rod_number == i_rod
            continue;
        end
        
        rr_j = voxel_list{rod_number};
        [rod_j,ind_j] = get_rod_from_centerline(rr_j,radius_threshold,stack);
        contact_indices = intersect(ind_i,ind_j);
                
        if numel(contact_indices) > 50 % ~isempty(contact_indices) % weak condition
            num_contact_voxels(end+1) = numel(contact_indices);
            highlight = [rod_i;rod_j];
            cuboid = get_bbox(size(stack),highlight);
            flash = get_line_img(cuboid([5,4,6])+[1,1,1],highlight - cuboid([2 1 3]) + [1 1 1]);
            
            cc = bwconncomp(flash);
            num_voxels_list = cellfun(@numel,cc.PixelIdxList);
            num_voxels_combined = (numel(ind_i) + numel(ind_j)) + numel(contact_indices)*2;
            
            
            contact_points = ind2sub2(size(stack),contact_indices);
            contacts_i(end+1) = rod_number;
            [dvec,d_min,sm_i,sm_j,e_i,e_j] = find_min_distance_btn_discrete_lines(rr_i,rr_j);
            
            p_1 = (sm_i + sm_j)/2;
            p_2 = mean(contact_points,1);
            n_ij = -dvec/norm(dvec);
            b_ij = cross(e_i,n_ij);
            t_ij = cross(n_ij,b_ij);
            Q_i0 = [n_ij,b_ij,t_ij];
            
            n_ji = -n_ij;
            b_ji = cross(e_j,n_ji);
            t_ji = cross(n_ji,b_ji);
            Q_j0 = [n_ji,b_ji,t_ji];
            
            i(end+1) = i_rod;
            j(end+1) = rod_number;
            p_i1(end+1,:) = p_1;
            p_i2(end+1,:) = p_2;
            s_i(end+1,:) = sm_i;
            s_j(end+1,:) = sm_j;
            Q_i(end+1,:) = Q_i0;
            Q_j(end+1,:) = Q_j0;
            
            if 0
                close all;
                tiledlayout(1,4)
                nexttile([1,3])
                plot3v(rod_i,'k.');hold on;
                plot3v(rod_j,'k.');axis equal;
                plot3v(contact_points,'ro');
                nexttile
                plot3v(contact_points,'ro');
                ;
            end
            
        end
    end
    num_contacts_i = numel(contacts_i);
    %     num_contact_voxels
    
    fprintf('%d: \t %d contact \t (%.2f sec elapsed).\n',i_rod,num_contacts_i,toc(t_start));
end

% close all;histogram(num_contact_voxels);

i = i';
j = j';

contact_table = table(i,j,p_i1,p_i2,s_i,s_j,Q_i,Q_j);

end