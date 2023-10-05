function splitted = bust_centerline(cl,threshold)
k = 1;
splitted = {};

t = bezier.eval(cl,1000);
dist_mat = pdist2(t,cl);

jj = min(dist_mat) < threshold;
ii = zeros(size(jj));
A = 1:numel(ii);

ii(strfind([0,jj(:)'],[0 1])) = 1;
idx = cumsum(ii).*jj;
contiguous_clusters = accumarray( idx(jj)',A(jj)',[],@(x){x'});

num_contiguous_clusters = numel(contiguous_clusters);

%     close all;
for i_cluster = 1:num_contiguous_clusters
    rr_here = cl(contiguous_clusters{i_cluster},:);
    splitted{end+1} = rr_here;
    %         plot3v(rr_here);hold on;axis equal;
    
end
k = k + 1;
if 0
    close all;
    plot3v(rr);hold on;axis equal
    plot3v(rr(min(dist_mat) < 1,:),'o');
end

end