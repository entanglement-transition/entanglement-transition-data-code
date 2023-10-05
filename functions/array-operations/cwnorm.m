function out = cwnorm(vec)    
    % column-wise norm
    out  = sqrt(sum(vec.^2,1));
end