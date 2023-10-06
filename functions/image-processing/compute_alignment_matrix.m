function [alignment_matrix,actual_distance_matrix] = compute_alignment_matrix(cl_list)
num_segments = numel(cl_list);
alignment_matrix = zeros(num_segments,num_segments);
actual_distance_matrix = zeros(num_segments,num_segments);
for ii = 1:num_segments
    seg_i = cl_list{ii};
    [~,ori_i,~] = get_line_coord(seg_i);
    for jj = ii+1:num_segments
        seg_j = cl_list{jj};
        [~,ori_j,~] = get_line_coord(seg_j);
        
        dist_mat = pdist2(seg_i,seg_j);
        [dvec,min_d,s1,s2,e1,e2]=find_min_distance_btn_discrete_lines(seg_i,seg_j);
        
        alignment_matrix(ii,jj) = (abs(dot(dvec/norm(dvec),ori_i)) + abs(dot(dvec/norm(dvec),ori_j)))/2;
        actual_distance_matrix(ii,jj) = min_d;
        
        %                 close all;plot3v(seg_i,'.-');hold on;plot3v(seg_j,'.-');axis equal;
        %                 plot3v(s1,'o','markersize',10);
        %                 plot3v(s2,'o','markersize',10);
    end
end
alignment_matrix = alignment_matrix + alignment_matrix';
actual_distance_matrix = actual_distance_matrix + actual_distance_matrix';

end