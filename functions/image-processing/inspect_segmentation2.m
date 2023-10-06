function I = inspect_segmentation2(stack,cline_list_,rad0)
%
stack_size = size(stack);
I_foreground = find(stack);
N_total = numel(I_foreground);
%%
k_stop = 100;
lwin = 10;
search_radius = 15;

% rad0 = 6;
hlen = 20;
marching_step = 5;

%%
cline_list = {};
tangent_list = {};
radius_list = [];
%%
k = 1;
tStart = tic;
N = numel(cline_list_);
all_found_points = cell(1,N);
ind_all = cell(1,N);
t_start = tic;
for i = 1:N % while removing stack points
    % get initial point
    %     [rad,centroid,orientation] = estimate_local_diameter(stack,rad0,hlen,search_radius);
    rr = round(cline_list_{i});
    [local_pc,ind] = cylinder_sweep_indices2(size(stack),rr,rad0);
    if isempty(ind)
        continue
    end
    I = stack(round(ind));
    ind = ind(I);
    ind_all{i} = ind;
    if mod(i,50) == 0
        fprintf('%d / %d \t %.2f sec\n',i,N,toc(t_start))
    end
end
I = vertcat(ind_all{:});
% 
% recon = zeros(size(stack),'logical');
% recon(I) = 1;
% 
% small = imresize3(stack,0.25);
% small_recon = imresize3(recon,0.25);


end