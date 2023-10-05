function new_cl_list = lengthen_centerlines_from_ends(zstack,cl_list,rad,cylinder_halflength)

marching_step = cylinder_halflength/2;
num_rods = numel(cl_list);
new_cl_list = cell(num_rods,1);
for i = 1:num_rods
    cl = cl_list{i};

    positive_axis = cl(end,:) - cl(1,:);
    positive_axis = positive_axis/norm(positive_axis);
    
    
    [found_point,found_axis,isend,local_volume,matching_error]...
        = locate_cylinder(zstack,round(cl(1,:)),rad,cylinder_halflength,positive_axis);
    [centerlines,tangents,to_add,matching_errors] =...
        segmentation_by_local_sampling_half(zstack,found_point,found_axis,rad,cylinder_halflength,marching_step);

    [found_point,found_axis,isend,local_volume,matching_error]...
        = locate_cylinder(zstack,round(cl(end,:)),rad,cylinder_halflength,-positive_axis);
    [centerlines2,tangents2,to_add2,matching_errors2] =...
        segmentation_by_local_sampling_half(zstack,found_point,found_axis,rad,cylinder_halflength,marching_step);

    new_cl_list{i} = reorder_centerline(vertcat(cl,centerlines,centerlines2));
%     close all;
%     plot3v(centerlines,'o-');hold on;
%     plot3v(centerlines2,'o-');
%     plot3v(cl,'k-');

    if mod(i,100) == 0
        fprintf('%d / %d completed.\n',i, num_rods);
    end
end

end


