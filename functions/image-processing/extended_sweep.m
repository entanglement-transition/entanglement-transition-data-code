function [sweep,indices] = extended_sweep(image_size,centerline,r)

% get voxel indices from inequalities
%     assert(mod(r,2) == 1)

if size(centerline,1) > 10
    first_five = centerline(1:10,:);
    [cen,ori,slist] = get_line_coord(first_five);
    
    first_extended = cen+ori.*(1:50)';
    
    close all;
    plot3v(centerline,'o');hold on;
    plot3v(first_extended,'ro-');
    ;
    
    last_five = centerline(end-4:end);
    
    
end


centerline = unique(round(centerline),'rows','stable');

ind = (1:2*r+1) - (r+1);
[X,Y,Z] = meshgrid(ind,ind,ind);
I = X.^2 + Y.^2 + Z.^2 <= r^2;
sphere_points = ind2sub2(size(I),find(I)) - (r+1);

num_centerline_points = size(centerline,1);
num_sphere_points = size(sphere_points,1);

all_points = zeros(num_centerline_points*num_sphere_points,3);
%     num_sphere_points*(1-1)+num_sphere_points+1;
for i = 1:num_centerline_points
    %         size( all_points(num_sphere_points*(i-1)+i:num_sphere_points*(i-1)+num_sphere_points,:) )
    %         size( sphere_points)
    
    all_points(num_sphere_points*(i-1)+1:num_sphere_points*(i-1)+num_sphere_points,:)...
        = centerline(i,:) + sphere_points;
end
sweep = unique(all_points,'rows','stable');
sweep = cut_to_matrix_size(stack_size,sweep);
indices = sub2ind2(stack_size,sweep);

end