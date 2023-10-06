function cl_list = remove_overlaps(cl_list,thick_rad)
nosed_list = segment_cells_to_nosed_list(cl_list);
all_points = nosed_list(:,1:3);
noses = nosed_list(:,4);

num_cls0 = numel(cl_list);

% overlap
overlap_segments = [];
for ii = 1:numel(cl_list)
    seg_i = cl_list{ii};
    fr_i = fit_rod(seg_i');
    
    [cen_i,ori_i,slist_i] = get_line_coord(seg_i);    
    slist = sum( ( all_points - cen_i).*ori_i,2 );
    dlist = rwnorm( all_points - ( cen_i + slist.*ori_i) );
    length_estimate = max(slist_i) - min(slist_i);
    indices_nearby_rods = dlist < thick_rad*2.5 & abs(slist) < length_estimate/2*1.5;
    nearby_rod_numberings = unique(noses(indices_nearby_rods));
    num_nearby_rods = numel(nearby_rod_numberings);


    for kk = 1:num_nearby_rods
        jj = nearby_rod_numberings(kk);
        if jj == ii
            continue;
        end

        seg_j = cl_list{jj};
        fr_j = fit_rod(seg_j');
        dist_mat = pdist2(fr_i.pts,fr_j.pts);
        if nnz(min(dist_mat) < thick_rad) > 750
            %                     close all;
            %                     plot3v(fr_i.pts,'.-');hold on;
            %                     plot3v(fr_j.pts,'.-');
            overlap_segments(end+1,:) = [jj,ii];
            ;
        elseif nnz(min(dist_mat') < thick_rad) > 750
            %                     close all;
            %                     plot3v(fr_i.pts,'.-');hold on;
            %                     plot3v(fr_j.pts,'.-');
            overlap_segments(end+1,:) = [ii,jj];
            ;
        end
    end
%     fprintf('%d / %d completed.\n',ii,numel(cl_list));
end

if isempty(overlap_segments)    
    return;
end
overlap_segments = unique(overlap_segments,'rows');
[C,ia,ic] = unique(sort(overlap_segments,2),'rows','stable');
unique_ones = overlap_segments(ia,1);

mutual_overlap = setdiff(1:size(overlap_segments,1),ia);
further_remove = [];
for i = mutual_overlap
    seg_i = cl_list{overlap_segments(i,1)};
    seg_j = cl_list{overlap_segments(i,2)};

    if sum(rwnorm(diff(seg_i))) > sum(rwnorm(diff(seg_j)))
        further_remove(end+1) = overlap_segments(i,2);
        unique_ones(unique_ones == overlap_segments(i,1)) = [];
    else
        further_remove(end+1) = overlap_segments(i,1);
        unique_ones(unique_ones == overlap_segments(i,2)) = [];
    end
end

to_remove = [unique_ones;further_remove'];

cl_list(to_remove) = [];
cl_list(cellfun(@(x) size(x,1),cl_list) < 2) = [];
% fprintf('%d no-overlap centerlines from %d original ones.\n',numel(cl_list),num_cls0 );

end


