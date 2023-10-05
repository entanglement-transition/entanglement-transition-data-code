function final_clusters = trim_centerlines(final_centerline_list,length_threshold,error_threshold)
% final_centerline_list = centerlines;
% length_threshold = 800;
% error_threshold = 1;

busting_threshold = 5;

[nn,ll,num_rods] = get_number_length(final_centerline_list);

close all;
histogram(ll);
print(gcf,'BeforeTrimming.png','-dpng','-r600');
%
error_list = zeros(num_rods,1);
for i = 1:num_rods
    cl = final_centerline_list{i};
    fr = fit_rod(cl');
    error_list(i) = mean(fr.err);
end

%
wrong_indices = error_list>error_threshold | ll' > length_threshold;
num_error_rods = nnz(wrong_indices);
%%
k = 1;
splitted = {};
for i = find(wrong_indices)'
    cl = final_centerline_list{i};
    splitted = join_two_cells(splitted,bust_centerline(cl,busting_threshold));
    k = k + 1;
    if 0
        close all;
        plot3v(rr);hold on;axis equal
        plot3v(rr(min(dist_mat) < 1,:),'o');
    end
end

%
busted = final_centerline_list;
busted(wrong_indices) = [];

busted = join_two_cells(busted,splitted);
num_rows = cellfun(@(x) size(x,1), busted);

%
[nn,ll,num_rods] = get_number_length(busted);
close all;
histogram(ll);
print(gcf,'Busted.png','-dpng','-r600');

disp('Busted the input centerlines!')
%
nosed_list = segment_cells_to_nosed_list(busted);
all_points = nosed_list(:,1:3);
noses = nosed_list(:,4);

%%
% error_threshold = 1;
friends = cell(num_rods,1);

feet_length_threshold = 10;
axial_length_threshold = 450;
%
tstart = tic;
for i_rod = 1:num_rods
    rr = busted{i_rod};
    %     len = sum(rwnorm(diff(rr)))
    
    [cen,ori,~] = get_line_coord(rr);
    slist = sum( ( all_points - cen).*ori,2 );
    dlist = rwnorm( all_points - ( cen + slist.*ori) );
    indices_nearby_rods = dlist < feet_length_threshold & abs(slist) < axial_length_threshold;
    
    nearby_rod_numberings = unique(noses(indices_nearby_rods));
    num_nearby_rods = numel(nearby_rod_numberings);
    local_nosed_list = nosed_list(indices_nearby_rods,:);
    
    %     close all;
    %     highlight = [];
    %     for i = 1:num_nearby_rods
    % %         plot3v(local_nosed_list( local_nosed_list(:,4) == nearby_rod_numberings(i), 1:3), '.-');
    % %         hold on;
    %         cl = local_nosed_list( local_nosed_list(:,4) == nearby_rod_numberings(i), 1:3);
    %         highlight = [highlight;cl];
    %     end
    %     visualize_highlighted_region(zstack,round(highlight));
    
    to_add = [];
    [~,ori0,~] = get_line_coord(rr);
    for i = 1:num_nearby_rods
        rod_number = nearby_rod_numberings(i);
        if rod_number == i_rod
            continue;
        end
        
        rr_join = join_two_rods(rr,busted{rod_number});
        fr = fit_rod(rr_join');
        
        [~,ori1,~] = get_line_coord(busted{rod_number});
        if mean(fr.err) < error_threshold & abs(dot(ori0,ori1)) > 0.95
            to_add(end+1) = i;
        end
    end
    friends{i_rod} = [i_rod,nearby_rod_numberings(to_add)'];
    
    if mod(i_rod,100) == 0
        fprintf('%d / %d friend groups are found. %.2f sec elapsed. \n',i_rod,num_rods,toc(tstart));
    end
    
    %     close all;
    %     figure('position',[0,0,1600,400])
    %     tiledlayout(1,2)
    %     nexttile
    %     for i = nearby_rod_numberings(to_add)'
    %         cl = local_nosed_list( local_nosed_list(:,4) == i, 1:3);
    %         plot3v(cl,'.-');hold on;
    %     end
    %     axis equal;
    %
    %     nexttile
    %     plot3v(busted([i_rod,nearby_rod_numberings(to_add)']),'.-');
    
end
%
friends_copy = friends;
num_friends = numel(friends_copy);
%%
groups = {};
while 1
    I = cellfun( @(x) any(ismember(friends_copy{1},x)),friends_copy);
    groups{end+1} = unique(horzcat(friends_copy{find(I)}));
    friends_copy(I) = [];
    numel(friends_copy);
    if isempty(friends_copy)
        break;
    end
end

%
num_groups = numel(groups);
ll = NaN(num_groups,1);
ll2 = NaN(num_groups,1);
for i = 1:num_groups
    rr = vertcat(busted{groups{i}});
    [~,~,slist] = get_line_coord(rr);
    [~,I] = sort(slist);
    rr = rr(I,:);
    
    fr = fit_rod(rr');
    if mean(fr.err) < 1
        ll(i) = sum(rwnorm(diff(rr)));
    else
        ll2(i) = sum(rwnorm(diff(rr)));
    end
end

%
clustered = groups(~isnan(ll));
not_clustered = groups(isnan(ll));
num_nans = nnz(isnan(ll));
not_clustered_numbers = find(isnan(ll));
optimal_clusters = {};

%%
small_length_threshold = 50;

t_start = tic;
cannot_solve = [];
for i = 1:num_nans    
    rod_numbers = groups{not_clustered_numbers(i)};
    cls_in_group = busted(rod_numbers);
    if isempty(cls_in_group)
        continue;
    end
    
    cls_in_group = remove_overlaps(cls_in_group,1.5);
    [num_points_list,individual_length_list,num_segments] = get_number_length(cls_in_group);    
    cls_in_group(individual_length_list < 2 | num_points_list < 3) = [];    
    if isempty(cls_in_group)
        continue;
    end    
    [num_points_list,individual_length_list,num_segments] = get_number_length(cls_in_group);
    
    % compute similarity/distance matrices
    alignment_matrix = compute_alignment_matrix(cls_in_group);
    error_matrix = compute_fitting_error_matrix(cls_in_group);
    dist_mat = 1./(alignment_matrix);
    
    K_max = min(10,numel(cls_in_group));
    K_min = 1;
    [bins,errors,C,sumd] = cluster_by_kmedoids(cls_in_group,dist_mat,individual_length_list,1,K_max);
    
    splitted = {};
    for ii = 1:max(bins)
        clustered_cl = vertcat(cls_in_group{bins == ii});
        clustered_cl = reorder_centerline(clustered_cl);
        splitted = join_two_cells(splitted,bust_centerline(clustered_cl,10));
    end    
    [~,ll,~] = get_number_length(splitted);
    splitted(ll < small_length_threshold) = [];
    
    %         close all;plot3v(splitted,'.');
    cls_in_group = splitted; % or just clustered?        
    if isempty(cls_in_group)
        continue;
    end
        
    [~,individual_length_list,~] = get_number_length(cls_in_group);
    num_segments = numel(cls_in_group);
    error_matrix = compute_fitting_error_matrix(cls_in_group);    
    [enemies_x,enemies_y] = ind2sub(size(error_matrix),find(error_matrix > 10));
    enemy_pairs = [enemies_x,enemies_y,enemies_y];
    
    % enemy_triple
    enemy_triples = [];
    if num_segments >= 3
        all_triples = nchoosek(1:num_segments,3);
        triple_error_list = zeros(size(all_triples,1),1);
        for ii = 1:size(all_triples,1)
            joined = vertcat( cls_in_group{all_triples(ii,:)} );
            [~,~,slist] = get_line_coord(joined);
            [~,I] = sort(slist);
            joined = joined(I,:);
            fr = fit_rod(joined');
            
            %         close all;
            %         plot3v(local_cell(all_triples(ii,:)),'o-');
            %         plot3v(fr.pts,'-');
            
            triple_error_list(ii) = mean(1/fr.r);
            
            if 0%mean(1/fr.r) > 1e-3
                %             fr.r
                close all;
                plot3v(cls_in_group(all_triples(ii,:)),'o-');
                plot3v(fr.pts,'-');
            end
        end
        enemy_triples = all_triples((triple_error_list > 1e-3),:);
    end    
    enemies = [enemy_pairs;enemy_triples];
    
    all_possible_partitions = partition(1:num_segments,4,enemies);    
    fprintf('Total %d ways of partitioning.\n',numel(all_possible_partitions));
    if numel(all_possible_partitions) > 1e5
        cannot_solve(end+1) = i;
        fprintf('Cannot solve this case. (%d)',i)
        continue;
        ;
    end
    
    k_max = numel(all_possible_partitions);
    error_list = zeros(k_max,1);
    combination_list = cell(k_max,1);
    min_length_list = zeros(k_max,1);;
    
    N = numel(all_possible_partitions);
    for ii = 1:N
        partition_ii = all_possible_partitions{ii};
        mean_error = 0;
        min_length = Inf;
        for j = 1:numel(partition_ii)
            
            cl_combined = vertcat(cls_in_group{partition_ii{j}});
            cl_combined = reorder_centerline(cl_combined);
                        
            fr = fit_rod(cl_combined');
            mean_error = mean_error + mean(fr.err); % consider using bezier curve            
            min_length = min(min_length,sum(individual_length_list(partition_ii{j})));
        end        
        error_list(k) = mean_error;
        min_length_list(k) = min_length;
        combination_list{k} = partition_ii;
        k = k + 1;
    end
    
    if k_max > 0
        [M,I] = min(error_list./min_length_list);
        optimal_partition = combination_list{I};
        N = numel(optimal_partition);
        
        %     close all;
        for ii = 1:N
            cl_combined = vertcat(cls_in_group{optimal_partition{ii}});
            cl_combined = reorder_centerline(cl_combined);

            
            optimal_clusters{end+1} = cl_combined;
            %         plot3v(rr_here,'.');hold on;axis equal;
        end
    elseif k_max == 0
        N = numel(cls_in_group);
        for ii = 1:N
            cl_combined = vertcat(cls_in_group{ii});
            cl_combined = reorder_centerline(cl_combined);

            optimal_clusters{end+1} = cl_combined;
            %         plot3v(rr_here,'.');hold on;axis equal;
        end
        optimal_partition = num2cell(1:num_segments);

    end
    
%     close all;
    for ii = 1:N
        cl_combined = vertcat(cls_in_group{optimal_partition{ii}});
        cl_combined = reorder_centerline(cl_combined);
%         plot3v(cl_combined,'.-'); hold on;
    end
%     ;axis equal
    ;
    if mod(i,10) == 0
        fprintf('%d collections resolved out of %d total collections in %.2f sec\n',i,num_nans,toc(t_start));
    end
    
end

%
num_clustered = numel(clustered);
clustered_by_fitting = cell(num_clustered,1);
for i = 1:num_clustered
    rr = vertcat(busted{clustered{i}});
    [~,~,slist] = get_line_coord(rr);
    [~,I] = sort(slist);
    rr = rr(I,:);
    clustered_by_fitting{i} = rr;
end
%
final_clusters = join_two_cells(clustered_by_fitting,optimal_clusters);
%
[~,ll,N] = get_number_length(final_clusters);

close all;
histogram(ll);
print(gcf,'AfterTrimming.png','-dpng','-r600');

end