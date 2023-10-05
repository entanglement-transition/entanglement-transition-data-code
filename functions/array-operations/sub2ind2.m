function ind = sub2ind2(sz,rr)
        ind = sub2ind(sz,rr(:,1),rr(:,2),rr(:,3));
%     if numel(sz) == 3
%         ind = sub2ind(sz,rr(:,1),rr(:,2),rr(:,3));
%     elseif numel(sz) == 2
%         ind = sub2ind(sz,rr(:,2),rr(:,1));
%     end
end