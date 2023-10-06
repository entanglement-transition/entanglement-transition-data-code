function [nb,previous_segments] = find_single_sided_neighbors(i,side_flag,segments_coordinates,labeled_list,num_smallest_segment,error_threshold,previous_segments)


N = numel(segments_coordinates);
rr = segments_coordinates{i};
R = num_smallest_segment * 2;
nb = [];

if side_flag == 1
    pivot_point = rr(1,:);
    other_end = rr(end,:);
elseif side_flag == 2
    pivot_point = rr(end,:);
    other_end = rr(1,:);
end

I = rwnorm(labeled_list(:,1:3) - pivot_point) < R;
labels1 = unique(labeled_list(I,4));
labels1 = setdiff(labels1,i);
error_list1 = zeros(1,numel(labels1));

for k = 1:numel(labels1)
    lb = labels1(k);
    rr2 = segments_coordinates{lb};
    [rr_new,e2e] = combine_segments(rr,rr2);    
    t = bezier.eval(rr_new,5*size(rr_new,1)); % param ***
% 
%     plot3v(e2e,'r-','linewidth',2);hold on;
%     plot3v(rr);
%     plot3v(rr2);
    
    errors = get_bezier_errors(rr,rr2,t);% for debug;
    error_list1(k) = errors.e5;
;
% close all;
% figure('Position',[-2400,0,500,500]);
%     plot3v(detached);
%     plot3v(rr,'r','linewidth',2);
%     hold on;
%     plot3v(rr2,'b','linewidth',2);
%     axis equal;
%     plot3(t(:, 1), t(:, 2), t(:,3),'-');%     plot3(t(TF, 1), t(TF, 2), t(TF,3),'o-');
%     text(rr2(1,1)+1,rr2(1,2)+1,rr2(1,3)+1,sprintf('%d, %.2f',lb,err));    
%     title(sprintf('%.2f',err*100));
    
    %
    %     for lb = [labels1(error_list1<5)']
    %         rr2 = segments_coordinates{lb};
    %         rr_new = combine_segments(rr,rr2);
    %
    %
    %         plot3v(rr,'r','linewidth',2);
    %         hold on;
    %         plot3v(rr2,'b','linewidth',2);
    %         axis equal;
    %
    %         t = bezier.eval(rr_new,5*size(rr_new,1));
    %         plot3(t(:, 1), t(:, 2), t(:,3),'.-');
    %
    %         dist_mat = pdist2(rr_new,t);
    %
    %         A = min(dist_mat,[],1);
    %         TF = isoutlier(A);
    %
    %         plot3v(rr2,'b','linewidth',2);
    %         text(rr2(1,1)+1,rr2(1,2)+1,rr2(1,3)+1,sprintf('%d_%.2f',lb,sum(A(TF))));
    %         plot3(t(:, 1), t(:, 2), t(:,3));
    %
    %     end
    %     axis equal;
    %     ;

end

[M,I] = sort(error_list1);
if numel(error_list1) > 4
    error_list1(I(5:end)) = [];
    labels1(I(5:end)) = [];
end

nb = labels1(error_list1<error_threshold)';

% [er,I] = sort(error_list1);
% I2 = error_list1<error_threshold;
% I = I & I2;
% nb = labels1( I(1:min(5,numel(labels1))) )';

new_nb = [];
num_labels = numel(nb);

if isempty(nb)
    nb = i;
    return;
else
    for k = 1:num_labels
        
        
        if numel(previous_segments) > 10
            continue;            
        end
        
        if ismember(lb,previous_segments)
            continue;
        end
        
        lb = nb(k);
        rr2 = segments_coordinates{lb};        
        d1 = norm(other_end-rr2(1,:));
        d2 = norm(other_end-rr2(end,:));
        d3 = norm(pivot_point-rr2(1,:));
        d4 = norm(pivot_point-rr2(end,:));
        
        [~,I] = min([d1,d2,d3,d4]);
                
        switch I
            case 1
                continue;
            case 2
                continue;
            case 3
                side_flag = 2;
                plot3v(pivot_point,'bo');
                plot3v(rr2,'.-');
                plot3v(rr2(1,:),'ko');
            case 4
                side_flag = 1;
                plot3v(pivot_point,'bo');
                plot3v(rr2,'.-');
                plot3v(rr2(end,:),'ko');
        end
        % when equality holds?
        
%         i
%         err        
%         rr2 = segments_coordinates{lb};
%         rr_new = combine_segments(rr,rr2);
%         t = bezier.eval(rr_new,5*size(rr_new,1)); % param ***
%         plot3v(t,'-','linewidth',2);
    
        previous_segments = unique([previous_segments,lb]);
        [nnb,ps] = find_single_sided_neighbors(lb,...
            side_flag,...
            segments_coordinates,...
            labeled_list,...
            num_smallest_segment,...
            error_threshold,...
            previous_segments);
        
        previous_segments = unique([previous_segments,ps]);
%         new_nb = [new_nb,nnb];
        new_nb = [new_nb,lb,nnb];
        
    end
end
nb = [nb,new_nb];
;

end