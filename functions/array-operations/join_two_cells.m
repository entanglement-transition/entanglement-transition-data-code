function c = join_two_cells(c1,c2)
assert(iscell(c1));
assert(iscell(c2));

N1 = numel(c1);
N2 = numel(c2);

N = N1 + N2;
c = cell(N,1);
for i = 1:N1
    c{i} = c1{i};
end

for i = 1:N2
    c{i + N1} = c2{i};
end

end