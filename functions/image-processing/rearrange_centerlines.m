function rr_out = rearrange_centerlines(rr,varargin)
if nargin > 1
    marching_step = varargin{1};
else
    marching_step = 10;
end
N = size(rr,1);
% [~,~,slist] = get_line_coord(rr);
% [sorted,i_sort] = sort(slist);
% rr_new = rr(i_sort,:);
rr_new = rr;
sorted = 1:N;

% [rr_new,ia,~] = unique(rr_new,'rows');

d = rwnorm(rr_new(2:end,:) - rr_new(1:end-1,:));
I = d < 1e-4;
rr_new(I,:) = [];
sorted(I)= [] ;
if size(rr_new,1) < 2
    rr_out = rr_new;
    return;
end

xx = sorted;
if norm(xx(end) - xx(1)) < 5
    rr_out = rr_new;
    return;
end
% xx = sorted(ia)';
% rr_out = interp1(xx,rr_new,linspace(xx(1),xx(end),60));
rr_out = interp1(xx,rr_new,xx(1):marching_step:xx(end));

end