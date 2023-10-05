
function fields = compute_contact_field(stack_size,centerlines,contact_locations,R,h)
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

system_height = max([labeled_centerline_edges(:,3);labeled_centerline_edges(:,6)]);
system_diameter = (max([labeled_centerline_edges(:,1);labeled_centerline_edges(:,4)]) + ...
    min([labeled_centerline_edges(:,2);labeled_centerline_edges(:,5)]) - ...
    min([labeled_centerline_edges(:,1);labeled_centerline_edges(:,4)]) - ...
    min([labeled_centerline_edges(:,2);labeled_centerline_edges(:,5)]) )/2;

system_volume = pi*( system_diameter )^2*( system_height )/4;
sphere_volume = 4*pi/3*R^3;


a = stack_size(1);
b = stack_size(2);
c = stack_size(3);
% a = round(system_diameter);
% b = round(system_diameter);
% c = round(system_height);

center_x = R:h:a-R;
center_y = R:h:b-R;
center_z = R:h:c-R;

num_x = size(center_x,2);
num_y = size(center_y,2);
num_z = size(center_z,2);

num_contacts = size(contact_locations,1);
expected_number_of_contacts = num_contacts*sphere_volume/system_volume;


fprintf('Size of the map: %d, %d, %d\n',num_x,num_y,num_z);
tic
tStart = tic;

contact_density_field = zeros(num_x,num_y,num_z);


for k = 1:num_z
    I_slab_ct = rwnorm(contact_locations(:,3) - center_z(k)) < 1.1*R;
    slab_ct = contact_locations(I_slab_ct,:);
    
    for i = 1:num_x
        for j = 1:num_y

            center = [center_x(i),center_y(j),center_z(k)];
            
            
            I_ct = rwnorm(slab_ct(:,1:3) - center) < R;
            num_contacts = nnz(I_ct);
            
%             contact_density_field(i,j,k) = num_contacts/expected_number_of_contacts;
            contact_density_field(i,j,k) = num_contacts;
            
        end
    end
    fprintf('Z-Layer: %d/%d \t Loop time: %.2f \t Elapsed time: %.2f\n',k,num_z,toc,toc(tStart));
    ;
end


fields.c = contact_density_field;

fields.cx = center_x;
fields.cy = center_y;
fields.cz = center_z;
fields.R = R;
fields.h = h;

end