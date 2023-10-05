function cuboid = get_bbox(image_size,rr)
%GET_BBOX Summary of this function goes here
%   Detailed explanation goes here
% TO DO: implement -0.5 +0.5 convention
x = rr(:,1);
y = rr(:,2);
z = rr(:,3);

lx = max(x)-min(x);
ly = max(y)-min(y);
lz = max(z)-min(z);

cuboid = [min(y) min(x) min(z) ly+1 lx+1 lz+1];

end

