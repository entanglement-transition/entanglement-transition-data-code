function cl_new = reorder_curly_centerlines(cl)

% identify endpoints
num_rows = size(cl,1);
dist_mat = pdist2(cl,cl);
diagonal = (1:num_rows) + (0:num_rows:(num_rows^2-num_rows));
dist_mat(diagonal) = Inf;

num_adjacent_points = zeros(num_rows,1);
for i = 1:num_rows
    num_adjacent_points(i) = nnz(dist_mat(i,:) < 1.8);
end

two_end_points = find(num_adjacent_points == 1);
assert(numel(two_end_points) == 2)

reordering = zeros(num_rows,1);

reordering(1) = two_end_points(1);
[~,I] = min(dist_mat(reordering(1),:));
reordering(2) = I;

for i = 2:num_rows-1
    i1 = reordering(i-1);
    i2 = reordering(i);

    [~,I] = sort(dist_mat(i2,:));

    i3 = setdiff(I(1:2),i1);

    d1 = dist_mat(i1,i2);
    d2 = dist_mat(i2,i3);
    d3 = dist_mat(i1,i3);

    assert(d3 <= d1 + d2)
    reordering(i+1) = i3;
end
assert(reordering(end) == two_end_points(2))

cl_new = cl(reordering,:);


end