function [rod,I_add_i] = get_rod_from_centerline(centerline,r,stack,varargin)
assert(~isempty(r));

[sweep_i,ind_i] = cylinder_sweep_indices2(size(stack),centerline,r);

% if size(r,1) == 1
%     [sweep_i,ind_i] = sphere_sweep_indices(size(stack),centerline,r);
% elseif size(r,1) > 1
%     [sweep_i,ind_i] = cylinder_sweep_indices2(size(stack),centerline,r);
% end
I_real_i = stack(ind_i);
I_add_i = ind_i(I_real_i);
rod = sweep_i(I_real_i,:);

% if nargin > 3
%     bbox = get_bbox(size(stack),rod);
%     crop = imcrop3(stack,bbox);
%     lifted = get_line_img(size(crop), rod - bbox([2 1 3]) + [1 1 1]);
%
%     % close all;volshow(lifted);
%     cc = bwconncomp(lifted);
%     num_pts = cellfun(@numel,cc.PixelIdxList);
%     [~,i_max] = max(num_pts);
%
%     lifted = zeros(size(lifted),'logical');
%     lifted(cc.PixelIdxList{i_max}) = 1;
%     skel = bwskel(lifted);
%
%     rod = ind2sub2(size(lifted),cc.PixelIdxList{i_max});
%     rod = rod + bbox([2 1 3]) - [1 1 1];
% end
% cc.PixelIdxList{i_max} = [];
% lifted( vertcat(cc.PixelIdxList{:}) ) = 0;
% lifted = imerode(lifted,strel('sphere',1));

end