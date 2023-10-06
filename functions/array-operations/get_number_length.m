function [nn,ll,N] = get_number_length(voxel_list)
    N = length(voxel_list);
    ll = zeros(1,N);
    nn = zeros(1,N);

    for i = 1:N
        rr = voxel_list{i};
%         fr = fit_rod(rr');
%         ll(i) = fr.len;
        [cen,ori,slist] = get_line_coord(rr);
        [~,I] = sort(slist);
        rr = rr(I,:);
        ll(i) = calculate_polygonal_line_length(rr);
        
        nn(i) = size(rr,1);
    end
    
end
