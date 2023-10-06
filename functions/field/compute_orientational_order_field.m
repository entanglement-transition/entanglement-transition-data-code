function fields = compute_orientational_order_field(stack_size,centerlines,R,h)
% R: radius of bounding spheres
% h: grid spacing

a = stack_size(1);
b = stack_size(2);
c = stack_size(3);

center_x = R:h:a-R;
center_y = R:h:b-R;
center_z = R:h:c-R;

num_x = size(center_x,2);
num_y = size(center_y,2);
num_z = size(center_z,2);

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

[XX,YY,ZZ] = meshgrid((1:2*R+1) - R,(1:2*R+1) - R,(1:2*R+1) - R);
I_sphere = XX.^2 + YY.^2 + ZZ.^2 < R^2;
num_sph_points = nnz(I_sphere);

fprintf('Size of the map: %d, %d, %d\n',num_x,num_y,num_z);
tic
tStart = tic;
number_field = zeros(num_x,num_y,num_z);

orientational_order_field = zeros(num_x,num_y,num_z);
number_field = zeros(num_x,num_y,num_z);
for k = 1:num_z
    I_slab_cl = rwnorm(labeled_centerline_edges(:,3) - center_z(k)) < 1.1*R;
    slab_cl = labeled_centerline_edges(I_slab_cl,:);
    
        
    for i = 1:num_x
        for j = 1:num_y
            center = [center_x(i),center_y(j),center_z(k)];
            I_cl = rwnorm(slab_cl(:,1:3) - center) < R;
            insider_labels = slab_cl(I_cl,7);
            insider_labels = unique(insider_labels);
            labeled_edges_bound = slab_cl(I_cl,:);            
            
            num_labels = numel(insider_labels);
%             linking_number_matrix = pdist2(labeled_edges_bound(:,1:6),labeled_edges_bound(:,1:6),@compute_lk);
% 
%             total_self_entanglement = 0;
%             for ii = 1:num_labels
%                 ind = labeled_edges_bound(:,7) == insider_labels(ii);
%                 e_self_matrix = pdist2(labeled_edges_bound(ind,1:6),labeled_edges_bound(ind,1:6),@compute_lk);
%                 e_self = sum(abs(e_self_matrix(:)),'omitnan');
%                 total_self_entanglement = total_self_entanglement + e_self;
%             end
%             
%             
            two_point_edges = zeros(num_labels,7);

            for i_lb = 1:num_labels
                e1 = labeled_edges_bound(labeled_edges_bound(:,7)==insider_labels(i_lb),1:3);
                e2 = labeled_edges_bound(labeled_edges_bound(:,7)==insider_labels(i_lb),4:6);
                p1 = e1(1,:);
                p2 = e2(end,:);                
                
                two_point_edges(i_lb,1:6) = [p1,p2];
                two_point_edges(i_lb,7) = insider_labels(i_lb);
            end
%             
% 
%             linking_number_matrix_2e = pdist2(two_point_edges(:,1:6),two_point_edges(:,1:6),@compute_lk);
%             total_self_entanglement_2e = 0;
%             for ii = 1:num_labels
%                 ind = two_point_edges(:,7) == insider_labels(ii);
%                 e_self_matrix = pdist2(two_point_edges(ind,1:6),two_point_edges(ind,1:6),@compute_lk);
%                 e_self = sum(abs(e_self_matrix(:)),'omitnan');
%                 total_self_entanglement_2e = total_self_entanglement_2e + e_self;
%             end
            
            number_field(i,j,k) = num_labels;
            
            tic
            orientational_order_field(i,j,k) = compute_orientational_order(two_point_edges(:,1:6));
            toc
                        
        end
    end
    fprintf('Z-Layer: %d/%d \t Loop time: %.2f \t Elapsed time: %.2f\n',k,num_z,toc,toc(tStart));
    ;
end


fields.n = number_field;
fields.s = orientational_order_field;

fields.cx = center_x;
fields.cy = center_y;
fields.cz = center_z;
fields.R = R;
fields.h = h;

end