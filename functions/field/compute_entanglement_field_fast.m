function fields = compute_entanglement_field_fast(zstack,centerlines,R,h)
% R: radius of bounding spheres
% h: grid spacing

[a,b,c] = size(zstack);

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

% num_contacts = size(contact_locations,1);
% system_height = max([labeled_centerline_edges(:,3);labeled_centerline_edges(:,6)]);
% system_diameter = (max([labeled_centerline_edges(:,1);labeled_centerline_edges(:,4)]) + ...
%     max([labeled_centerline_edges(:,2);labeled_centerline_edges(:,5)]) - ...
%     min([labeled_centerline_edges(:,1);labeled_centerline_edges(:,4)]) - ...
%     min([labeled_centerline_edges(:,2);labeled_centerline_edges(:,5)]) )/2;
% 
% system_volume = pi*( system_diameter )^2*( system_height )/4;
% sphere_volume = 4*pi/3*R^3;
% expected_number_of_contacts = num_contacts*sphere_volume/system_volume;

[XX,YY,ZZ] = meshgrid((1:2*R+1) - R,(1:2*R+1) - R,(1:2*R+1) - R);
I_sphere = XX.^2 + YY.^2 + ZZ.^2 < R^2;
num_sph_points = nnz(I_sphere);

fprintf('Size of the map: %d, %d, %d\n',num_x,num_y,num_z);
tic
tStart = tic;
number_field = zeros(num_x,num_y,num_z);

% contact_density_field = zeros(num_x,num_y,num_z);
% volume_fraction_field = zeros(num_x,num_y,num_z);
entanglement_field = zeros(num_x,num_y,num_z);
% entanglement_field_2e = zeros(num_x,num_y,num_z);
% orientational_order_field = zeros(num_x,num_y,num_z);

for k = 1:num_z
    I_slab_cl = rwnorm(labeled_centerline_edges(:,3) - center_z(k)) < 1.1*R;
    slab_cl = labeled_centerline_edges(I_slab_cl,:);
    
%     I_slab_ct = rwnorm(contact_locations(:,3) - center_z(k)) < 1.1*R;
%     slab_ct = contact_locations(I_slab_ct,:);
    
    for i = 1:num_x
        for j = 1:num_y
            center = [center_x(i),center_y(j),center_z(k)];
            I_cl = rwnorm(slab_cl(:,1:3) - center) < R;

            insider_labels = slab_cl(I_cl,7);
            insider_labels = unique(insider_labels);
            labeled_edges_bound = slab_cl(I_cl,:);
            
%             I_ct = rwnorm(slab_ct(:,1:3) - center) < R;
%             num_contacts = nnz(I_ct);
            
            leftmost = max(min(size(zstack),[center([2,1,3])-R]),[0.5,0.5,0.5]);
            cuboid = [leftmost,2*R*ones(1,3)];
            crop = imcrop3(zstack,cuboid);
            
            num_labels = numel(insider_labels);
            
            linking_number_matrix = pdist2(labeled_edges_bound(:,1:6),labeled_edges_bound(:,1:6),@compute_lk);
            total_self_entanglement = 0;
            for ii = 1:num_labels
                ind = labeled_edges_bound(:,7) == insider_labels(ii);
                e_self_matrix = pdist2(labeled_edges_bound(ind,1:6),labeled_edges_bound(ind,1:6),@compute_lk);
                e_self = sum(abs(e_self_matrix(:)),'omitnan');
                total_self_entanglement = total_self_entanglement + e_self;
            end
            
            
               
%             number_field(i,j,k) = num_labels;
%             volume_fraction_field(i,j,k) = nnz(crop & I_sphere)/num_sph_points;
%             orientational_order_field(i,j,k) = compute_orientational_order(labeled_edges_bound(:,1:6));
%             contact_density_field(i,j,k) = num_contacts/expected_number_of_contacts;
%             entanglement_field(i,j,k) = sum(abs(linking_number_matrix(:))) - total_self_entanglement;
            entanglement_field_2e(i,j,k) = sum(abs(linking_number_matrix_2e(:))) - total_self_entanglement_2e;
            
        end
    end
    fprintf('Z-Layer: %d/%d \t Loop time: %.2f \t Elapsed time: %.2f\n',k,num_z,toc,toc(tStart));
    ;
end


% fields.n = number_field;
% fields.phi = volume_fraction_field;
% fields.s = orientational_order_field;
% fields.c = contact_density_field;
% fields.e = entanglement_field;
fields.e2 = entanglement_field_2e;

fields.cx = center_x;
fields.cy = center_y;
fields.cz = center_z;
fields.R = R;
fields.h = h;

end