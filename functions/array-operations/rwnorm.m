function out = rwnorm(x)
% row-wise norm
out = sqrt(sum(x.^2,2));
end