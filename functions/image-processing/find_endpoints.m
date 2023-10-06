function endpoints = find_endpoints(rr)

    distance_matrix = pdist2(rr,rr);
    
    N = size(rr,1);
    
    endpoints = [];
    for i = 1:N
        [sorted,I] = sort(distance_matrix(i,:));
        
        if sorted(2) ~= sorted(3)
            endpoints(end+1) = i;
        end 

    end
    ;

end