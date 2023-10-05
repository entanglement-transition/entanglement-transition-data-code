function [centerlines,tangents,found_points,matching_errors] = segmentation_by_local_sampling_half(stack,initial_point,initial_axis,rad,lwin,marching_step)
% marching_step = lwin/4; % TUNE HERE

% first head "left"
current_point = initial_point;
current_axis = -initial_axis;

centerlines = [];
tangents = [];
found_points = {};
matching_errors = [];

k = 1;
while 1 % going "left"
    [centerline_point,tangent_vector,isend,local_volume,matching_error] = locate_cylinder(stack,current_point,rad,lwin,current_axis);
    if isend
        centerlines(end+1,:) = centerline_point;
        tangents(end+1,:) = tangent_vector;
        found_points{end+1} = local_volume;
        matching_errors(end+1) = matching_error;
        break;
    end
    centerlines(end+1,:) = centerline_point;
    tangents(end+1,:) = tangent_vector;
    found_points{end+1} = local_volume;
    matching_errors(end+1) = matching_error;
    
    % update
    current_point = round(centerline_point + (marching_step)*tangent_vector);
    current_axis = tangent_vector;
    
    k = k + 1;
    if k > 500
        break;
    end
end
current_point = initial_point;
current_axis = initial_axis; % head right

end

