function connected_segments = follow_branches(i,R,segments_coordinates,labeled_list)

rr = segments_coordinates{i};

I = rwnorm(labeled_list(:,1:3) - rr(1,:)) < R;
labels1 = unique(labeled_list(I,4));

connected_segments = [];    
if numel(labels1) == 1
    return;
else
    for lb = labels1
        connected_segments = [connected_segments,follow_branches(i,R,segments_coordinates,labeled_list)];
    end
end


end