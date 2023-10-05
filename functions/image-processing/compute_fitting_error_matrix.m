function error_matrix = compute_fitting_error_matrix(cl_list)
num_segments = numel(cl_list);
error_matrix = zeros(num_segments,num_segments);
for ii = 1:num_segments
    seg_i = cl_list{ii};
    for jj = ii+1:num_segments
        seg_j = cl_list{jj};
        
        joined = join_two_rods(seg_i,seg_j);
        [~,~,slist] = get_line_coord(joined);
        [~,I] = sort(slist);
        joined = joined(I,:);
        
        fr = fit_rod(joined');
        error_matrix(ii,jj) = mean(fr.err);
    end
end
error_matrix = error_matrix + error_matrix';
end