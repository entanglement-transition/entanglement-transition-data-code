function [bins,errors,C,sumd] = cluster_by_kmedoids(local_cell,dist_mat,length_list,K_min,K_max)
num_segments = numel(local_cell);
errors = NaN(K_max,1);
%
for K = K_min:K_max
    [bins,C,sumd] = kmedoids((1:num_segments)',K,'Distance',@(I,J) dist_mat(I,J));
    
    mean_errors = zeros(max(bins),1);
    length_sum = zeros(max(bins),1);    
    for ii = 1:max(bins)
        cl_combined = vertcat(local_cell{bins == ii});
        cl_combined = reorder_centerline(cl_combined);
        fr = fit_rod(cl_combined');
        
        mean_errors(ii) = mean(fr.err);
        length_sum(ii) = sum(length_list(bins == ii));        
    end
    title(sprintf('%.8f',sum(mean_errors)/min(length_sum)));
    
    errors(K) = sum(mean_errors.^2)/mean(length_sum.^2);
    
end
%
[~,K_opt] = min(errors,[],'omitnan');
[bins,C,sumd] = kmedoids((1:num_segments)',K_opt,'Distance',@(I,J) dist_mat(I,J));
;


end