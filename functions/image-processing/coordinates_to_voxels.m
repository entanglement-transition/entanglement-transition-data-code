function [voxels,corner] = coordinates_to_voxels(rr)

x = rr(:,1);
y = rr(:,2);
z = rr(:,3);

lx = max(x)-min(x);
ly = max(y)-min(y);
lz = max(z)-min(z);

corner = [min(x),min(y),min(z)];%-[0.5,0.5,0.5];
image_size = [lx+1,ly+1,lz+1];

rr = rr - corner + [1,1,1];
% rr = cut_to_matrix_size(image_size,rr);

ind = sub2ind2(image_size,rr);

voxels = zeros(image_size,'logical');
voxels(ind) = 1;


end