function [cen,ori,slist] = get_line_coord(rr)
    cen = mean(rr,1);
    [~,~,V] = svd(rr-cen);
    ori = V(:,1)';
    slist = sum((rr-cen).*ori,2);
end