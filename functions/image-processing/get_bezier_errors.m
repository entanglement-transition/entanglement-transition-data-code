function errors = get_bezier_errors(rr,rr2,t)

[rr_new,e2e] = combine_segments(rr,rr2);

dist_mat0 = pdist2(e2e,t);
B0 = min(dist_mat0,[],2);


dist_mat = pdist2(rr_new,t);
A = min(dist_mat,[],1);
B = min(dist_mat,[],2);

% close all;
% plot3v(rr_new,'b-'); hold on;
% plot3v(rr_new(B > 1,:),'ro');
% plot3v(t,'k-');
% plot3v(t(A>5,:),'c--');
 
TF_outlier = isoutlier(A);
TF_large = A > 1;

detached = t(TF_large,:);
jj = A > 1;
ii = zeros(size(jj));
B = 1:numel(ii);

ii(strfind([0,jj(:)'],[0 1])) = 1;
idx = cumsum(ii).*jj;
contiguous_clusters = accumarray( idx(jj)',B(jj)',[],@(x){x'});

if ~isempty(contiguous_clusters)
    cluster_sizes = cellfun(@(x) numel(x), contiguous_clusters);
    [~,I] = max(cluster_sizes);

    detached = t(contiguous_clusters{I},:);
    diff1 = detached(3:end,:) - detached(2:end-1,:);
    diff2 = detached(2:end-1,:) - detached(1:end-2,:);
    curvature_sum = mean(rwnorm(diff2 - diff1)) * 100;
else
    curvature_sum = 0;
end

errors.e1 = sum(A(TF_outlier));
errors.e2 = sum(A(TF_large));

errors.e3 = curvature_sum;
errors.e4 = sum(B(B>1));
errors.e5 = sum(B0);
errors.e6 = sum(A(TF_outlier));
errors.e7 = sum(A(TF_outlier));


end