
function I_cylinder = cylinder_matching_indices(rr,x,rad,hlen)
% x y z azi ele
rr0 = [x(1) x(2) x(3)];
[u,v,w] = sph2cart(x(4),x(5),1);
ax = [u,v,w];
slist = sum( (rr-rr0).*ax, 2);
dlist = rwnorm( rr - (rr0 + slist.*ax ) );

I_cylinder  = (dlist < rad) & (abs(slist) < hlen);

end