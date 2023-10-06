function out = reorder_centerline(cl)

[~,~,slist] = get_line_coord(cl);
[~,I] = sort(slist);
out = cl(I,:);

end