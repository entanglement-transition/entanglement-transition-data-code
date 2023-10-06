function nb = find_neighbor_tree(i,segments_coordinates,labeled_list,num_smallest_segment,error_threshold)


N = numel(segments_coordinates);
rr = segments_coordinates{i};
R = num_smallest_segment * 2;
I = rwnorm(labeled_list(:,1:3) - rr(1,:)) < R;
labels1 = unique(labeled_list(I,4));
labels1 = setdiff(labels1,i);

I = rwnorm(labeled_list(:,1:3) - rr(end,:)) < R;
labels2 = unique(labeled_list(I,4));
labels2 = setdiff(labels2,i);

error_list1 = zeros(1,numel(labels1));
for k = 1:numel(labels1)
    lb = labels1(k);
    rr2 = segments_coordinates{lb};
    rr_new = combine_segments(rr,rr2);
    t = bezier.eval(rr_new,5*size(rr_new,1));
    
    dist_mat = pdist2(rr_new,t);
    A = min(dist_mat,[],2);
    TF = isoutlier(A);
    
    error_list1(k) = sum(A(TF));
    
    
end

error_list2 = zeros(1,numel(labels2));
for k = 1:numel(labels2)
    lb = labels2(k);
    rr2 = segments_coordinates{lb};
    rr_new = combine_segments(rr,rr2);
    t = bezier.eval(rr_new,5*size(rr_new,1));
    
    dist_mat = pdist2(rr_new,t);
    A = min(dist_mat,[],2);
    TF = isoutlier(A);
    
    error_list2(k) = sum(A(TF));
end

nb = [labels1(error_list1<error_threshold)',labels2(error_list2<error_threshold)'];

if isempty(nb)
    ;
    return
else
    for lb = nb
        new_nb = find_neighbor_tree(lb,segments_coordinates,labeled_list,num_smallest_segment,error_threshold);
        nb = [nb,new_nb];
    end
end


% visual
% close all;
% figure('Position',[-2400,0,500,500]);
% plot3v(rr,'ro','linewidth',2);
% hold on;
%
% for lb = [labels1(error_list1<5)',labels2(error_list2<5)'] % ***
%     rr2 = segments_coordinates{lb};
%     rr_new = combine_segments(rr,rr2);
%
%
%     plot3v(rr,'r','linewidth',2);
%     hold on;
%     plot3v(rr2,'b','linewidth',2);
%     axis equal;
%
%     t = bezier.eval(rr_new,5*size(rr_new,1));
%     plot3(t(:, 1), t(:, 2), t(:,3),'.-');
%
%     dist_mat = pdist2(rr_new,t);
%
%     A = min(dist_mat,[],2);
%     TF = isoutlier(A);
%
%     plot3v(rr2,'b','linewidth',2);
%     text(rr2(1,1)+1,rr2(1,2)+1,rr2(1,3)+1,sprintf('%d_%.2f',lb,sum(A(TF))));
%     plot3(t(:, 1), t(:, 2), t(:,3));
%
% end
% axis equal;


end