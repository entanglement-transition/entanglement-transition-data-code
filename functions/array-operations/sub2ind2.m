function ind = sub2ind2(sz,rr)

rr = round(rr);
rr = max(rr,[1,1,1]);
rr = min(rr,sz);

        ind = sub2ind(sz,rr(:,1),rr(:,2),rr(:,3));
%     if numel(sz) == 3
%         ind = sub2ind(sz,rr(:,1),rr(:,2),rr(:,3));
%     elseif numel(sz) == 2
%         ind = sub2ind(sz,rr(:,2),rr(:,1));
%     end
end