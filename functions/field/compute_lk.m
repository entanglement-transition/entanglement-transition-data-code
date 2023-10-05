function out = compute_lk(ZI,ZJ)
% Panagiotou and Kauffman 2020 (and its ref [27]: Banchoff)
num_row = size(ZJ,1);
p1 = repmat(ZI(:,1:3),[num_row, 1]);
p2 = repmat(ZI(:,4:6),[num_row, 1]);
q1 = ZJ(:,1:3);
q2 = ZJ(:,4:6);

r11 = p1 - q1;
r12 = p1 - q2;
r21 = p2 - q1;
r22 = p2 - q2;

n1 = get_normal(r11,r12);
n2 = get_normal(r12,r22);
n3 = get_normal(r22,r21);
n4 = get_normal(r21,r11);

sgn = sign(dot(cross(p2-p1,q2-q1,2),r11,2));
t1 = clip(dot(n1,n2,2),-1,1);
t2 = clip(dot(n2,n3,2),-1,1);
t3 = clip(dot(n3,n4,2),-1,1);
t4 = clip(dot(n4,n1,2),-1,1);

mag = asin( t1 ) + asin( t2 ) + asin( t3 ) + asin( t4 );
out = sum(sgn.*mag, 2) / (4*pi);
end