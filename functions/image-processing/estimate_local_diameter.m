function [rad,centroid,orientation,local_volume,cyl_ind,score,num_skips]...
    = estimate_local_diameter(stack,I_search,found_ind,rad0,hlen,search_radius)
% rad0: initial value for radius
% hlen: half length of cylinder
% search radius: radius of local search volume
stack_size = size(stack);
num_skips = 0;

cline_list = {};
tangent_list = {};

% search_radius = 15;
rad = rad0;
mag = @(x) sqrt(sum(x.^2,2));

local_volume = [];
centroid = [];
orientation = [];
cyl_ind = [];
score = 0;
% I_search = [];

num_skips = 1;
tStart = tic;

while 1 % while finding initial point and axis
    %         ind = randsample(I_search,1);
    if num_skips > 100 % parameter
        return;
    end
        
    pick = randsample(numel(I_search),1);    
    ind = I_search(pick);
    if ~stack(ind)        
        break;
    end
    
    current_point = ind2sub2(stack_size,ind);
    
    current_axis = [1,0,0];
    [sph,ind] = sphere_sweep_indices(stack_size,current_point,search_radius); % also need to tune this number; especially when we deal with larger diameter
    search = sph(stack(ind),:);
    lb = point_segmentation(search);
        
    is_good = 0;
    scores = size(1,max(lb));
    
    for l = 1:max(lb)        
        local_volume = search(lb==l,:);
        if size(local_volume,1) < 0.4 * pi*rad0.^2*search_radius
            scores(l) = 0;
            num_skips = num_skips + 1;
            continue;
        end
        
        stats = get_principal_axis_length(local_volume);
        V = stats.EigenVectors;
        
        hlen = max(stats.PrincipalAxisLength)/sqrt(2);
        current_point = mean(local_volume,1);
        current_axis = V(:,1)';
        
        [u,v,w] = cart2sph(current_axis(1),current_axis(2),current_axis(3));
        x0 = [current_point u v];
        cyl_ind = cylinder_matching_indices(local_volume,x0,rad0+1,hlen);
        score = nnz(cyl_ind) /size(local_volume,1);
        scores(l) = score;        
    end
    [max_score,where] = max(scores);
    
    if max_score == 0
        num_skips = num_skips + 1;
        continue;
    end
    
    
    local_volume = search(lb == where,:);
    stats = get_principal_axis_length(local_volume);
    V = stats.EigenVectors;
    hlen = max(stats.PrincipalAxisLength)/sqrt(2);
    current_point = mean(local_volume,1);
    current_axis = V(:,1)';
    
    [u,v,w] = cart2sph(current_axis(1),current_axis(2),current_axis(3));
    x0 = [current_point u v];
    cyl_ind = cylinder_matching_indices(local_volume,x0,rad0+1,hlen);
    score = nnz(cyl_ind) /size(local_volume,1);
    
    if ismember(local_volume,found_ind)
        num_skips = num_skips + 1;
        continue;
    end
    
    if 0%size(local_volume,1) < 1000
        close all;
        plot3v(local_volume,'.'); hold on;
        plot3v([current_point;current_point + hlen*current_axis],'r-');
        
        title(size(local_volume,1));axis equal
        ;
    end
    
    %         [u,v,w] = cart2sph(current_axis(1),current_axis(2),current_axis(3));
    %         x0 = [current_point u v];
    
    opts = optimset('Display','none','MaxFunEvals',5000);
    num_trials = 20;
    rad_trials = linspace(1,rad0,num_trials);
    evals = zeros(1,num_trials);
    
    for tt = 1:num_trials
        rad = rad_trials(tt);
        [x,fval,exitflag,output] = fminsearch(@(x) cylinder_matching(local_volume,x,rad,hlen),x0,opts);
        evals(tt) = -fval;
    end
    
    [~,i_max] = max(evals);
    rad = rad_trials(i_max);
    cyl_ind = cylinder_matching_indices(local_volume,x,rad,hlen);
    score = nnz(cyl_ind)/size(local_volume,1);
    [u,v,w] = sph2cart(x(4),x(5),1);
    tangent_vector = [u,v,w];
    cen = x(1:3);
    tvec = sum((local_volume-cen).*tangent_vector,2)*tangent_vector;
    tmag = mag(tvec);
    rvec = (local_volume-cen) - tvec;
    rmag = mag(rvec);
    rg = sqrt(sum(rmag.^2)/numel(rmag) )*sqrt(2);
    %         score = nnz(ind_in)/numel(cyl_ind);
    num_skips = num_skips + 1;    
    
    if 0
        close all;
        yj.plot3(local_volume,'.');hold on;
        yj.plot3(current_point,'r.');
        yj.plot3([current_point;current_point+ V(:,1)'*hlen],'linewidth',2);axis equal;
    end
    
    if score > 0.75 % parameter
        
        %             nnz(ind_in)/numel(cyl_ind);
        %             figure
        %             plot(tmag,rmag,'.');
        %             hold on;
        %             plot(tmag(cyl_ind,:),rmag(cyl_ind,:),'.');
        %             plot([tmag(1),tmag(end)],[rg,rg],'-','linewidth',2);
        
        if 0 % DEBUG
            [X,Y,Z] = cylinder(rad);
            Qf = rotMat([0,0,1]',tangent_vector');
            Z = 2*(Z-0.5)*hlen;
            temp=[X(:),Y(:),Z(:)]*Qf';
            sz=size(X);
            X=reshape(temp(:,1),sz) + x(1);
            Y=reshape(temp(:,2),sz) + x(2);
            Z=reshape(temp(:,3),sz) + x(3);
            
            figure;
            plot3v(local_volume,'.');
            title(sprintf('score: %.2f \t radius: %.2f',score,rad))
            hold on;
            plot3v(local_volume(cyl_ind,:),'r.');
            plot3v(current_point,'r.');
            plot3v(x(1:3),'ro');
            plot3v([x(1:3);x(1:3)+tangent_vector*hlen],'linewidth',2);
            axis equal
            surf(X,Y,Z,'facecolor','w','facealpha',0.2)
            axis equal
            %                 close all;
            %                 ind_in = rmag(cyl_ind) < rg;
            
        end
        
        orientation = tangent_vector;
        centroid = round([x(1),x(2),x(3)]);
        
        break;        
    end
    
    
    num_skips = num_skips + 1;
    
end

end

%%
function out = cylinder_matching(rr,x,rad,hlen)
% x y z azi ele
rr0 = [x(1) x(2) x(3)];
[u,v,w] = sph2cart(x(4),x(5),1);
ax = [u,v,w];
slist = sum( (rr-rr0).*ax, 2);
dlist = rwnorm( rr - (rr0 + slist.*ax ) );

I_cylinder  = (dlist < rad);
% I_cylinder  = (dlist < rad) & (abs(slist) < hlen);
% out = -nnz(I_cylinder);
out = pi*rad.^2*hlen - nnz(I_cylinder);

end

function I_cylinder = cylinder_matching_indices(rr,x,rad,hlen)
% x y z azi ele
rr0 = [x(1) x(2) x(3)];
[u,v,w] = sph2cart(x(4),x(5),1);
ax = [u,v,w];
slist = sum( (rr-rr0).*ax, 2);
dlist = rwnorm( rr - (rr0 + slist.*ax ) );

I_cylinder  = (dlist < rad) & (abs(slist) < hlen);

end

function rot = rotMat(b,a,alpha)
% ROTMAT returns a rotation matrix that rotates unit vector b to a
%
%   rot = rotMat(b) returns a d x d rotation matrix that rotate
%   unit vector b to the north pole (0,0,...,0,1)
%
%   rot = rotMat(b,a ) returns a d x d rotation matrix that rotate
%   unit vector b to a
%
%   rot = rotMat(b,a,alpha) returns a d x d rotation matrix that rotate
%   unit vector b towards a by alpha (in radian)
%
%    See also .

% Last updated Nov 7, 2009
% Sungkyu Jung


[s1 s2]=size(b);
d = max(s1,s2);
b= b/norm(b);
if min(s1,s2) ~= 1 || nargin==0 , help rotMat, return, end

if s1<=s2;    b = b'; end

if nargin == 1;
    a = [zeros(d-1,1); 1];
    alpha = acos(a'*b);
end

if nargin == 2;
    alpha = acos(a'*b);
end
if abs(a'*b - 1) < 1e-15; rot = eye(d); return, end
if abs(a'*b + 1) < 1e-15; rot = -eye(d); return, end

c = b - a * (a'*b); c = c / norm(c);
A = a*c' - c*a' ;

rot = eye(d) + sin(alpha)*A + (cos(alpha) - 1)*(a*a' +c*c');
end
