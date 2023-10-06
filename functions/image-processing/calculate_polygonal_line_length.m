function len = calculate_polygonal_line_length(rr)
    len = sum(sqrt(sum((rr(2:end,:) - rr(1:end-1,:)).^2,2)),1);
end