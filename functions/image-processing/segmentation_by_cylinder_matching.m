function [segmented_centerlines,segmented_tangents,segmented_radii,segmented_length] = segmentation_by_cylinder_matching(zstack,radius_0)

cylinder_halflength = radius_0*3; %20
search_radius = cylinder_halflength; %10
marching_step = ceil(search_radius/2); %5

%%
stack_size = size(zstack);
I_foreground = find(zstack);
N_total = numel(I_foreground);
N_search = size(I_foreground,1);
N_remain = N_search;

%%
segmented_centerlines = {};
segmented_tangents = {};
segmented_radii = [];
segmented_length = [];
found_ind = [];
N_found = 0;
k = 0;
l = 0;
tStart = tic;
while 1 % while removing stack points
    % get initial point
    tstart = tic;
    [rad,centroid,orientation,local_volume,cyl_ind,score,num_skips] = estimate_local_diameter(zstack,I_foreground,found_ind,radius_0,cylinder_halflength,search_radius);
    num_skips;
    tm = toc(tstart);
    if isempty(orientation)
        continue;
        ;
    end
    if tm>5
        fprintf('Took some time to segment: %.6f sec\n',tm);
    end

    % marching
    [centerlines,tangents,to_add,matching_errors] = segmentation_by_local_sampling(zstack,centroid,orientation,rad,cylinder_halflength,marching_step);
    cl = rearrange_centerlines(centerlines,cylinder_halflength/2);
    [sweep,sweep_ind] = cylinder_sweep_indices2(size(zstack),cl,rad*1.05);
    true_sweep = zstack(sweep_ind);

    if any(true_sweep)
        k = k + 1;
        zstack(sweep_ind(true_sweep)) = 0;

        N_found = N_found + nnz(true_sweep);
        N_remain = N_total - N_found;

        segmented_centerlines{end+1} = centerlines;
        segmented_tangents{end+1} = tangents;
        segmented_radii(end+1) = rad;
        segmented_length(end+1) = sum(rwnorm(diff(centerlines)));

        fprintf('%d, \t radius = %.2f, \t score = %.2f, \t %.2f sec elapsed \n',k,rad,score,toc(tStart));
        fprintf('length of just segmented rod: %.2f \t average length: %.2f \n',segmented_length(end),mean(segmented_length));
        fprintf('Remaining: %d, \t Removed: %d, \t Initial: %d, \t Complete: %2.4f %% \n',N_remain,N_found,N_total,100*(N_found)/N_total);

    end
    l = l + 1;
    if mod(l,10) == 0
        volume_threshold = pi*rad^2*(cylinder_halflength*1);
        zstack = clear_zstack(zstack,volume_threshold);

        N_remain = nnz(zstack);
        N_found = N_total - N_remain;
    end

    

    if N_remain/N_total < 0.10
%     if sum(segmented_length) > 625 * 16
        break;
    end
end

end