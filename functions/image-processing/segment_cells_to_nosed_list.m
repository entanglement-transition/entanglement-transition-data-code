function nosed_list = segment_cells_to_nosed_list(segment_cell)
    num_segments = length(segment_cell);
    all_points = vertcat(segment_cell{:});
    num_segments = length(segment_cell);
    num_all_points = size(all_points,1);
    nose = zeros(num_all_points,1);
    initial_cells = cell(num_segments,1);
    l0 = 1;    
    for i = 1:num_segments
        l = size(segment_cell{i},1);
        nose( l0 : l0 + l - 1) = i;
        l0 = l0 + l;
    end
    nosed_list = [all_points nose];
end