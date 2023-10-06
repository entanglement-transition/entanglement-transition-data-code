function [sweep,indices] = cylinder_sweep_indices2(stack_size,cen,r)
    
    % get voxel indices from inequalities
%     assert(mod(r,2) == 1)
    cen = unique(round(cen),'rows','stable');
    if size(cen,1) < 2
        sweep = [];
        indices = [];
        return;
    end
    num_centerline_points = size(cen,1);
%     num_cylinder_points = size(cylinder_points,1);    
%     all_points = zeros(num_centerline_points*num_cylinder_points,3);
%     num_sphere_points*(1-1)+num_sphere_points+1;
    all_points_cell = cell(1,num_centerline_points);
    for i = 1:num_centerline_points-1        
        h = norm(cen(i+1,:) - cen(i,:)); % parameter
        
        meshgrid_size = max(h,r);
        
        ind = (1:2*meshgrid_size+1) - (meshgrid_size+1);
        [X,Y,Z] = meshgrid(ind,ind,ind);

        tt = (cen(i+1,:) - cen(i,:))/norm(cen(i+1,:) - cen(i,:));
        
        M1 = X.^2 + Y.^2 + Z.^2;
        M2 = (tt(1)*X + tt(2)*Y + tt(3)*Z).^2;
        I = M1 - M2 <= r^2;
        cylinder_points = ind2sub2(size(I),find(I)) - (h+1);
        if i == num_centerline_points-1
            
            ;
        end
%         cylinder_points = ind2sub2(size(I),find(I));        
%         all_points(num_cylinder_points*(i-1)+1:num_cylinder_points*(i-1)+num_cylinder_points,:)...
%             = cen(i,:) + cylinder_points;
        all_points_cell{i} = cen(i,:) + round(cylinder_points);
    end
    all_points = vertcat(all_points_cell{:});
    sweep = unique(all_points,'rows','stable');    
    sweep = cut_to_matrix_size(stack_size,sweep);
    indices = sub2ind2(stack_size,sweep);
    
%     close all;plot3v(round(cylinder_points) + cen(i,:),'.');axis equal;hold on;plot3v(cen(i,:),'o','markersize',10);plot3v(cen(end,:),'o');plot3v([cen(i+1,:);cen(i,:)]);plot3v(all_points);axis equal;
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

