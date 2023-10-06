function [cl_list,good_segments,scores,error_list] = segmentation_by_bwskel(zstack,radius_0,score_threshold,derror_threshold)


skel = bwskel(zstack);
bran = bwmorph3(skel,'branchpoints');

% se = strel('sphere',2);
% bran = imdilate(bran,se);
% lines = bwmorph3(skel&~bran,'clean'); % remove isolated voxels

lines = skel&~bran;
cc = bwconncomp(lines);

ind_list = cc.PixelIdxList;
num_elements = cellfun(@numel, ind_list);
ind_list(num_elements<20) = [];

num_cc = numel(ind_list);
error_list = zeros(num_cc,1);
cl_list = cell(num_cc,1);

scores = NaN(num_cc,1);
tstart = tic;
for i = 1:num_cc
    ind = ind_list{i};
    cl = ind2sub2(size(zstack),ind);

    [cen,ori,slist] = get_line_coord(cl);
    [~,I] = sort(slist);
    cl = cl(I,:);
    cl_list{i} = cl;

    dlist = rwnorm(cl - (cen + ori.*slist));
    error_list(i) = mean(dlist);

    if mean(dlist) < 100 & numel(ind) > 5;
        % TO DO: smooth cl otherwise its sweep will be very erratic
        [sweep,sweep_ind] = cylinder_sweep_indices2(size(zstack),cl,radius_0);
        true_sweep = zstack(sweep_ind);
        if all(~true_sweep)
            continue
            ;
        end

        rod = sweep(true_sweep,:);
        stats = get_principal_axis_length(rod);
        V = stats.EigenVectors;
        hlen = max(stats.PrincipalAxisLength)/sqrt(2);
        current_point = mean(rod,1);
        current_axis = V(:,1)';

        [u,v,w] = cart2sph(current_axis(1),current_axis(2),current_axis(3));
        x0 = [current_point u v];
        cyl_ind = cylinder_matching_indices(rod,x0,radius_0+1,hlen);
        score = nnz(cyl_ind) /size(rod,1);
        scores(i) = score;
    end
    
    if mod(i,100) == 0
        fprintf('%d / %d processed. %.2f sec elapsed.\n',i,num_cc,toc(tstart));
    end
end
good_segments = scores > score_threshold & error_list < derror_threshold;

end