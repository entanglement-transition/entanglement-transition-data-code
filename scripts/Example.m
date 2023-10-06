restoredefaultpath
setup

mkdir('../results')
%%
container_radius = 600;
container_height = 600;
container_volume = pi*container_radius^2*container_height;

edge_length = 600;
rod_diameter = 6;

Z = 3;
num_rods = Z*container_volume/rod_diameter/edge_length^2; % num_rods/container_volume*rod_diameter*edge_length^2 = z

number_density = num_rods/container_volume;
number_density*edge_length^2*rod_diameter

random_edges = generate_intersecting_random_rods_in_cylinder(num_rods,container_radius,container_height,edge_length);

%%
rod_radius = rod_diameter/2;
N = size(random_edges,1);
distance_matrix = Inf(N,N);
for i = 1:N
    edge_i = random_edges(i,:);
    distance_lower_bound = extended_line_distances(edge_i,random_edges);
    
    j_select = find(distance_lower_bound < 10*rod_radius)';
    j_select = setdiff(j_select,i);
    
    for j = j_select
        edge_j = random_edges(j,:);
        distance_matrix(i,j) = distance_between_edges(edge_i,edge_j);
    end
end

%%

[contact_info,distance_matrix] = get_contact_info_for_edges(random_edges,rod_radius);
[~,ia,~] = unique(sort([contact_info.i,contact_info.j],2),'row');
contacts = contact_info(ia,:);
size(contacts)

g = graph(contacts.i,contacts.j);
deg = degree(g);

%%

N = size(random_edges,1);
rod_image = zeros(2*container_radius,2*container_radius,container_height+1,'logical');
for i = 1:N
    r1 = random_edges(i,1:3) + [container_radius+1,container_radius+1,1];
    r2 = random_edges(i,4:6) + [container_radius+1,container_radius+1,1];
    
    valid_voxels = unique(round(linspacev(r1,r2,1000)),'row');
    
    I = sub2ind2(size(rod_image),valid_voxels);
    rod_image(I) = 1;
end

%% generate 3D array of voxel iamage

zstack = imdilate(rod_image,strel('sphere',rod_radius));

[cl_list,good_segments] = segmentation_by_bwskel(zstack,ceil(rod_radius),0.8,0.5);
trimmed = trim_centerlines(cl_list,800,1);
final_centerline_list = remove_overlaps(trimmed,rod_radius);
for i = 1:numel(final_centerline_list)
    rr = final_centerline_list{i};
    smoothen = smoothdata(rr);
    
    final_centerline_list{i} = rr;
    
%     close all;
%     plot3v(smoothen,'r-');
%     hold on;
%     plot3v(rr(1:10:end,:),'b.-');
    ;
end

radius_0 = rod_radius;
cylinder_halflength = radius_0 * 4;
final_centerline_list = lengthen_centerlines_from_ends(zstack,final_centerline_list,radius_0,cylinder_halflength);
trimmed = trim_centerlines(final_centerline_list,750,1);
fitted_centerlines = cell(size(trimmed));
for i = 1:numel(trimmed)
    rr = trimmed{i};
    fr = fit_rod(rr');
        
    fitted_centerlines{i} = fr.pts;
end
toc

%%
shifted_random_edges = random_edges + [container_radius+1,container_radius+1,1,container_radius+1,container_radius+1,1];    
N = size(shifted_random_edges,1);

segmented_edges = get_edges(fitted_centerlines);

centerlines = cell(N,1);
score_list = zeros(N,1);

tic
for i = 1:N
    edge_i = shifted_random_edges(i,:);
    rr1 = [edge_i(1:3);edge_i(4:6)];    
    
    distance_lower_bound = extended_line_distances(edge_i,segmented_edges);
    
    j_select = find(distance_lower_bound < 10*rod_radius)';
    j_select = setdiff(j_select,i);
    
    min_distances = zeros(size(j_select));
    sum_of_pairwise_distances = zeros(size(j_select));
    k = 1;
    for j = j_select
        
        rr2 = fitted_centerlines{j};        
        [dvec,min_d,s1,s2,e1,e2] = find_min_distance_btn_discrete_lines(rr1,rr2);        
        min_distances(k) = min_d;        
        
        dist_mat = pdist2(rr1,rr2);
        
        sum_of_pairwise_distances(k) = sum(min(dist_mat,[],2));
        
        k = k + 1;
    end
    [M,I] = min(sum_of_pairwise_distances);
    score_list(i) = M;
    centerlines{i} = fitted_centerlines{j_select(I)};
    toc
end

%%
set_figure(6,5);
histogram(score_list(score_list<2),'normalization','probability');
hold on;
% histogram(score_list(score_list<2),'normalization','probability');
xlabel('$\bar{\epsilon}_\mathrm{centerline}$ (pixels)')
ylabel('Probability')
print(gcf,'../results/centerline_error.png','-dpng','-r600');

%%
tic
contacts = get_contact_info(centerlines,zstack,rod_radius+0.5);
toc

%%
fitted_centerlines = cell(size(centerlines));
for i = 1:N
    rr = centerlines{i};
    fr = fit_rod(rr');
        
    fitted_centerlines{i} = fr.pts;
end

%%
score_list = zeros(N,1);
for i = 1:N
    rr1 = [shifted_random_edges(i,1:3);shifted_random_edges(i,4:6)];
    rr2 = fitted_centerlines{i};
    
    M = size(rr2,1);
    min_distances = zeros(M,1);
    for j = 1:M
        min_distances(j) = find_min_distance_btn_point_line(rr2(j,:),rr1);        
    end
    score_list(i) = mean(min_distances);
end
%%
[d,dist_vec,contact_site] = distance_between_edges(edge_i,edge_j);

%%
set_figure(6,5)
histogram(score_list,linspace(0,1,30),'normalization','probability')
xlabel('$\bar{\epsilon}_\mathrm{centerline}$ (pixels)')
ylabel('Probability');
print(gcf,'../results/segmentation_error.png','-dpng','-r600');


set_figure(6,5);
histogram(score_list,linspace(1,100,30),'normalization','probability')
xlabel('$\bar{\epsilon}_\mathrm{centerline}$ (pixels)')
ylabel('Probability');
print(gcf,'../results/segmentation_error2.png','-dpng','-r600');

%%
generated_edges = shifted_random_edges;
generated_image = zstack;

% found_contacts = contacts;
found_centerlines = centerlines;

%%
[generated_contacts,~]= get_contact_info_for_edges(generated_edges,rod_radius);

%%
generated_nodes = [generated_contacts.i,generated_contacts.j];
srtd = sort(generated_nodes,2);
[~,ia] = unique(srtd,'row');
reduced_generated_contacts = generated_contacts(ia,:);

%
found_nodes = [found_contacts.i,found_contacts.j];
srtd = sort(found_nodes,2);
[~,ia] = unique(srtd,'row');
reduced_found_contacts = found_contacts(ia,:);
%%
matched_indices = [];
num_unmatched = 0;
error_in_contact_position = [];

for i = reduced_generated_contacts.i'
    j_list = reduced_generated_contacts(reduced_generated_contacts.i == i,:).j';
    
    for j = j_list        
        matched = reduced_found_contacts(reduced_found_contacts.i == i & reduced_found_contacts.j == j,:);
        
        if isempty(matched)
            num_unmatched = num_unmatched + 1;
            continue
        else
            matched_indices(end+1,:) = [i,j];
            p1 = reduced_found_contacts(reduced_found_contacts.i == i & reduced_found_contacts.j == j,:).p_i1;
            p2 = reduced_generated_contacts(reduced_generated_contacts.i == i & reduced_generated_contacts.j == j,:).p;
            
            error_in_contact_position(end+1) = norm(p1-p2);            
        end        
    end
end
%%
set_figure(6,5);
histogram(error_in_contact_position(error_in_contact_position<1),25,'normalization','probability');
xlabel('$\epsilon_\mathrm{contact}$ (pixels)')
ylabel('Probability')
print(gcf,'../results/contact_error.png','-dpng','-r600');

%%
function [image_stacks,left_corner,right_corner] = edges_to_volume(edges,rod_radius)

N = numel(edges_to_volume,1);

left_corners = zeros(N,3);
right_corners = zeros(N,3);

voxellist_to_export = cell(N,1);

fat_rods = cellfun(@(x) sphere_sweep(x,rod_radius*2),voxel_list,'uniformoutput',false);

voxellist_to_export{i} = fat_rods;

voxels = round(vertcat(fat_rods{:}));
left_corners(i,:) = min(voxels);
right_corners(i,:) = max(voxels);

left_corner = min(left_corners);
right_corner = max(right_corners);

left_corner = left_corner - [10,10,10];
right_corner = right_corner + [20,20,20];
image_size = right_corner - left_corner;

image_stacks(N).stack = zeros(image_size,'logical');


voxel_cell = cellfun(@(x) unique(round(x)-left_corner,'rows','stable'),voxellist_to_export{i},'uniformoutput',false);
voxel_array = vertcat(voxel_cell{:});
stack = zeros(image_size,'logical');
ind = sub2ind2(image_size,voxel_array);
stack(ind) = 1;
image_stacks(i).stack = stack;



end

function [contact_info,distance_matrix]= get_contact_info_for_edges(random_edges,rod_radius)
% tolerance?

i = [];
j = [];
p = [];


N = size(random_edges,1);
distance_matrix = Inf(N,N);
for i_rod = 1:N
    edge_i = random_edges(i_rod,:);
    distance_lower_bound = extended_line_distances(edge_i,random_edges);
    
    j_select = find(distance_lower_bound < 10*rod_radius)';
    %     neighbor_edges = random_nonintersecting_rod_edges(distance_lower_bound < 10*rod_radius,:);
    
    for j_rod = j_select
        edge_j = random_edges(j_rod,:);
        [d,dvec,contact_point] = distance_between_edges(edge_i,edge_j);
        
        if d < rod_radius*2
            
            distance_matrix(i_rod,j_rod) = d;
            
            i(end+1) = i_rod;
            j(end+1) = j_rod;
            p(end+1,:) = contact_point;
            
        end
    end
end

i = i'; j = j';
contact_info = table(i,j,p);

% nnz(distance_matrix < rod_radius*2)/size(random_edges,1)/2
% contact_info.average_contact_number = nnz(distance_matrix < rod_radius*2)/size(random_edges,1)/2;
% contact_info.


end