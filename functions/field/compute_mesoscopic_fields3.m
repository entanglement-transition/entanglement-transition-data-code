function fields = compute_mesoscopic_fields3(centerlines,R,h,rod_radius)
% R: radius of bounding spheres
% h: grid spacing

rods_in_cylinder = vertcat(centerlines{:});
max_values = max(rods_in_cylinder);
min_values = min(rods_in_cylinder);

stack_size = ceil((max_values - min_values)/100)*100;
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

system_height = max([labeled_centerline_edges(:,3);labeled_centerline_edges(:,6)]);
system_diameter = (max([labeled_centerline_edges(:,1);labeled_centerline_edges(:,4)]) + ...
    min([labeled_centerline_edges(:,2);labeled_centerline_edges(:,5)]) - ...
    min([labeled_centerline_edges(:,1);labeled_centerline_edges(:,4)]) - ...
    min([labeled_centerline_edges(:,2);labeled_centerline_edges(:,5)]) )/2;

system_volume = pi*( system_diameter )^2*( system_height )/4;
sphere_volume = 4*pi/3*R^3;

[XX,YY,ZZ] = meshgrid((1:2*R+1) - R,(1:2*R+1) - R,(1:2*R+1) - R);
I_sphere = XX.^2 + YY.^2 + ZZ.^2 < R^2;
num_sph_points = nnz(I_sphere);

fprintf('Size of the map: %d, %d, %d\n',num_x,num_y,num_z);
tic
tStart = tic;
number_field = zeros(num_x,num_y,num_z);

volume_fraction_field = zeros(num_x,num_y,num_z);
entanglement_field = zeros(num_x,num_y,num_z);
entanglement_field_2e = zeros(num_x,num_y,num_z);
orientational_order_field = zeros(num_x,num_y,num_z);

for k = 1:num_z
    I_slab_cl = rwnorm(labeled_centerline_edges(:,3) - center_z(k)) < 1.1*R;
    slab_cl = labeled_centerline_edges(I_slab_cl,:);

    for i = 1:num_x
        for j = 1:num_y
            center = [center_x(i),center_y(j),center_z(k)];
            I_cl = rwnorm(slab_cl(:,1:3) - center) < R;

            if any(I_cl)

                insider_labels = slab_cl(I_cl,7);
                insider_labels = unique(insider_labels);
                labeled_edges_bound = slab_cl(I_cl,:);
                
                num_labels = numel(insider_labels);
                linking_number_matrix = pdist2(labeled_edges_bound(:,1:6),labeled_edges_bound(:,1:6),@compute_lk);

                total_self_entanglement = 0;
                for ii = 1:num_labels
                    ind = labeled_edges_bound(:,7) == insider_labels(ii);
%                     e_self_matrix = pdist2(labeled_edges_bound(ind,1:6),labeled_edges_bound(ind,1:6),@compute_lk);
                    e_self_matrix = 1;
                    e_self = sum(abs(e_self_matrix(:)),'omitnan');
                    total_self_entanglement = total_self_entanglement + e_self;
                end

                two_point_edges = zeros(num_labels,7);
                for i_lb = 1:num_labels
                    e1 = labeled_edges_bound(insider_labels==insider_labels(i_lb),1:3);
                    e2 = labeled_edges_bound(insider_labels==insider_labels(i_lb),4:6);
                    p1 = e1(1,:);
                    p2 = e2(end,:);

                    two_point_edges(i_lb,1:6) = [p1,p2];
                    two_point_edges(i_lb,7) = insider_labels(i_lb);
                end

                linking_number_matrix_2e = pdist2(two_point_edges(:,1:6),two_point_edges(:,1:6),@compute_lk);
                total_self_entanglement_2e = 0;
                for ii = 1:num_labels
                    ind = two_point_edges(:,7) == insider_labels(ii);
%                     e_self_matrix = pdist2(two_point_edges(ind,1:6),two_point_edges(ind,1:6),@compute_lk);
                    e_self_matrix = 1;
                    e_self = sum(abs(e_self_matrix(:)),'omitnan');
                    total_self_entanglement_2e = total_self_entanglement_2e + e_self;
                end

                lengths = rwnorm(two_point_edges(:,1:3) - two_point_edges(:,4:6));
                phi = sum(pi*rod_radius^2*lengths)/(4*pi/3*R^3);

                number_field(i,j,k) = num_labels;
                volume_fraction_field(i,j,k) = phi;
                orientational_order_field(i,j,k) = compute_orientational_order2(two_point_edges(:,1:6));
                entanglement_field(i,j,k) = sum(abs(linking_number_matrix(:))) - total_self_entanglement;
                entanglement_field_2e(i,j,k) = sum(abs(linking_number_matrix_2e(:))) - total_self_entanglement_2e;
            end

        end
    end
    fprintf('Z-Layer: %d/%d \t Loop time: %.2f \t Elapsed time: %.2f\n',k,num_z,toc,toc(tStart));
    ;
end


fields.n = number_field;
fields.phi = volume_fraction_field;
fields.s = orientational_order_field;
fields.e = entanglement_field;
fields.e2 = entanglement_field_2e;

fields.cx = center_x;
fields.cy = center_y;
fields.cz = center_z;
fields.R = R;
fields.h = h;

end