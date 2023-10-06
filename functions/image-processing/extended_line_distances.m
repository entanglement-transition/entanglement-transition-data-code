function d = extended_line_distances(edge,e2)

if isempty(e2)
    d = Inf;
    return
end
N = size(e2,1);
e1 = repmat(edge,[N,1]);

p1 = e1(:,1:3);
q1 = e1(:,4:6);
p2 = e2(:,1:3);
q2 = e2(:,4:6);

d1 = q1 - p1;
d2 = q2 - p2;

r = p1 - p2;

a = dot(d1, d1,2);
b = dot(d1, d2,2);
c = dot(d2, d2,2);
e = dot(d1, r,2);
f = dot(d2, r,2);

det = a.*c - b.^2;

I = abs(det) < 1e-6;
% I = (1 - a.*c/b.^2) < 1e-6;

s = (b.*f - c.*e)./det;
t = (a.*f - b.*e)./det;

dist_vec = r + s.*d1 - t.*d2;
d = rwnorm(dist_vec);
d(I) = 0;

end