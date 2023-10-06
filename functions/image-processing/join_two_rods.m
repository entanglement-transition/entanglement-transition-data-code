function rr = join_two_rods(r1,r2)
    rr = [r1;r2];
    [~,~,slist] = get_line_coord(rr);
    
    [~,I] = sort(slist);
    rr = rr(I,:);
end