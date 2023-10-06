function [found_point,found_axis,isend,local_volume,matching_error] = locate_cylinder(stack,query_point,rad,hlen,query_axis)
search_radius = hlen;

isend = 0;
current_point = query_point;
current_axis = query_axis;
local_volume = [];
matching_error = Inf;

if 0 % figure script
    % ==== %
    % stack_points = ind2sub2(size(stack),find(stack));
    % fh = figure(99);
    % plot3v(stack_points(1:20:end,:),'.');hold on
    % ==== %
end


% cut local volume to cylinder shape
cline = current_point + (-hlen:hlen)'.*current_axis;
[sph,ind] = sphere_sweep_indices(size(stack),cline,search_radius); % also need to tune this number; especially when we deal with larger diameter

search = sph(stack(ind),:);
lb = point_segmentation(search);
which_label = lb(find(ismember(search,round(current_point),'rows')));
if isempty( which_label ) % exit through wall can happen...
    isend = 1;
    found_point = current_point;
    found_axis = current_axis;
    matching_error = Inf;
    return;
end
local_volume = search(lb == which_label,:);

% update
%     previous_point = current_point;
previous_axis = current_axis;

[u,v,w] = cart2sph(current_axis(1),current_axis(2),current_axis(3));
x0 = [current_point u v];
x = fminsearch(@(x) cylinder_matching(local_volume,x,rad,hlen),x0);
[matching_error,ind_in] = cylinder_matching(local_volume,x,rad,hlen);
local_volume = local_volume(ind_in,:);

if matching_error > 0.5
    isend = 1;
    found_point = current_point;
    found_axis = current_axis;
    return;
    ;
end

current_point = x(1:3);
[u,v,w] = sph2cart(x(4),x(5),1);
current_axis = [u,v,w];
current_axis = current_axis*sign( dot(previous_axis,current_axis) );

% test if this is end point
test = round(current_point + (0:1.5*hlen)'.*current_axis); % is 10 a magic number? - should tune at some point
test = cut_to_matrix_size(size(stack),test);
test_ind = sub2ind2(size(stack), test);

if ~all(stack(test_ind))
    current_point = test(max(find(stack(test_ind))),:);
    isend = 1;
    
%     close all;    
%     plot3v(local_volume,'.');hold on;
%     plot3v(current_point,'o','linewidth',2,'markersize',10);
%     plot3v([current_point;current_point+10*current_axis]);
%     title(sprintf('%.4f',matching_error)); axis equal;
    ;
end

found_point = current_point;
found_axis = current_axis;

% if backward
% if dot( (current_point - query_point), query_axis) < 0
% if norm(current_point - query_point) < hlen/10
%     isend = 1;
%     found_point = query_point;
%     found_axis = query_axis;
%     return;
% end

% cylinder_radius = 3.5;
% cylinder_halflength = lwin;
% [pp,slist,dlist] = project_points_to_line(local_volume,current_point,current_axis);
% I_cylinder  = (dlist < cylinder_radius) & (abs(slist) < cylinder_halflength);
% num_cyl_pts = nnz(I_cylinder);

if 0 % in case of fat cylinder
    % find cylinder by iteration
    stats.PrincipalAxisLength(2)
    
    %         cylinder_radius = 3.5;
    %         cylinder_halflength = 20;
    %         [pp,slist,dlist] = project_points_to_line(local_volume,current_point,local_axis);
    %         I_cylinder  = (dlist < cylinder_radius) & (abs(slist) < cylinder_halflength);
    %
    %         close all
    %         plot3v(search,'ro');hold on
    %         plot3v(local_volume,'o');axis equal;
    %         plot3v(current_point + (-200:200)'.*local_axis,'-','linewidth',2);
    %         plot3v(local_volume(I_cylinder,:),'o');
    
    %         [u,v,w] = sph2cart(local_axis(1),local_axis(2),local_axis(3));
    [u,v,w] = cart2sph(current_axis(1),current_axis(2),current_axis(3));
    
    x0 = [current_point u v];
    x = fminsearch(@(x) cylinder_matching(local_volume,x),x0);
    %         options = optimset('Display','iter','PlotFcns',@optimplotfval);
    %         x = fminsearch(@(x) cylinder_matching(local_volume,x),x0,options);
    current_point = x(1:3);
    [u,v,w] = sph2cart(x(4),x(5),1);
    local_axis = [u,v,w];
    
    cylinder_radius = 3.5;
    cylinder_halflength = hlen;
    [pp,slist,dlist] = project_points_to_line(local_volume,current_point,local_axis);
    I_cylinder  = (dlist < cylinder_radius) & (abs(slist) < cylinder_halflength);
    
    close all;
    figure('position',[800 100 1500 1000])
    plot3v(local_volume,'o');hold on;axis equal
    plot3v(current_point,'o','markersize',15,'linewidth',3);
    plot3v(local_volume(I_cylinder,:),'ro')
    plot3v(mean(local_volume,1),'o','markersize',15,'linewidth',3);
    plot3v(current_point + (-20:20)'.*local_axis);
    plot3v(pp,'r-');
    %     title(sprintf('%d',cylinder_correlation(i)))
    ;
    
    cylinder_stats = get_principal_axis_length(local_volume(I_cylinder,:));
    cylinder_axis = cylinder_stats.EigenVectors(:,1)';
    error_angle = norm(cross(local_axis,cylinder_axis));
    % find cylinder by matching
    % position & angle both!
    % or use current position and change angle only
end

% case handling: if distance difference is too large
if 0%error_angle > 30
    angle_diff0 = norm(cross(current_axis,stat0.EigenVectors(:,1)'))
    angle_diff = norm(cross(stats.EigenVectors(:,1)',stat0.EigenVectors(:,1)'))
    
    if angle_diff > 0
        set(fh,'position',[2300 100 1500 1000])
        %         figure('position',[800 100 1500 1000])
        plot3v(local_volume,'o','linewidth',2);hold on;axis equal
        plot3v(local_volume(I_cylinder,:),'o','linewidth',2);hold on;axis equal
        plot3v(search,'k.')
        
        plot3v(current_point,'bo','markersize',15,'linewidth',3);
        plot3v(mean(local_volume,1),'ro','markersize',15,'linewidth',3);
        
        plot3v(current_point + (-20:20)'.*stat0.EigenVectors(:,1)');
        plot3v(current_point + (-20:20)'.*current_axis);
        
        plot3v(mean(local_volume,1) + (-20:20)'.*stats.EigenVectors(:,1)');
        title(sprintf('dist diff: %.4f, cyl corr: %d',dist_diff,nnz(I_cylinder)));
        
        ;
    end
    ;
end

% update
%     previous_point = current_point;
if 0 % figure script
    % =================================================================== %
    %     N = 100; % TUNE HERE
    %     orientation_list = randvonMisesFisherm(3, N, 10, stats.EigenVectors(:,1)')';
    %     cylinder_radius = 4;
    %     cylinder_halflength = lwin;
    %     for i = 1:N
    %         [pp,slist,dlist] = project_points_to_line(local_volume,current_point,orientation_list(i,:));
    %         I_cylinder  = (dlist < cylinder_radius) & (abs(slist) < cylinder_halflength);
    %         cylinder_correlation(i) = nnz(I_cylinder); % TUNE HERE
    %
    %         stat0 = get_principal_axis_length(local_volume(I_cylinder,:));
    %
    %         norm(cross(orientation_list(i,:),stat0.EigenVectors(:,1)'))
    %
    %         close all;
    %         figure('position',[800 100 1500 1000])
    %         plot3v(local_volume,'o');hold on;axis equal
    %         plot3v(current_point,'o','markersize',15,'linewidth',3);
    %         plot3v(local_volume(I_cylinder,:),'ro')
    %         plot3v(mean(local_volume,1),'o','markersize',15,'linewidth',3);
    %         plot3v(current_point + (-20:20)'.*stat0.EigenVectors(:,1)');
    %         plot3v(pp,'r-');
    %         title(sprintf('%d',cylinder_correlation(i)))
    %         ;
    %     end
    %     [M,I] = max(cylinder_correlation);
    %     max_correlation = M;
    %     max_orientation = orientation_list(I,:);
    %
    %     if M < pi*cylinder_radius^2*(2*cylinder_halflength)*0.8;
    %         ;
    %     end
    
    % =================================================================== %
end

% 
% if (1-dot(previous_axis, current_axis)^2 < 0.01)...
%         & (num_cyl_pts < 150000)...
%         %1 - 0.98^2 = 0.04 % TUNE HERE
%     
%     found_point = current_point;
%     found_axis = current_axis;
%     
%     break;
% end
% 


end

function [out,I_cylinder] = cylinder_matching(rr,x,rad,hlen)
% x y z azi ele
rr0 = [x(1) x(2) x(3)];
[u,v,w] = sph2cart(x(4),x(5),1);
ax = [u,v,w];
slist = sum( (rr-rr0).*ax, 2);
dlist = rwnorm( rr - (rr0 + slist.*ax ) );

% I_cylinder  = (dlist < rad) & (abs(slist) < hlen);
% out = -nnz(I_cylinder);

% % I_cylinder  = (dlist < rad);
I_cylinder  = (dlist < rad) & (abs(slist) < hlen);
% out = -nnz(I_cylinder);
out = (pi*rad.^2*(2*hlen) - nnz(I_cylinder))/(pi*rad.^2*(2*hlen));

end