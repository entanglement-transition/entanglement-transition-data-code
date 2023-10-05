function h = plot3v2(rr,varargin)
    if ismatrix(rr) & ~iscell(rr)
        assert(size(rr,2) == 3);
        h = plot3(rr(:,1),rr(:,2),rr(:,3),varargin{:});    
    elseif iscell(rr) & ~isempty(rr)        
        N = numel(rr);       
        for i = 1:N
            plot3(rr{i}(:,1),rr{i}(:,2),rr{i}(:,3),varargin{:});
            hold on
            text(rr{i}(1,1)+1,rr{i}(1,2),rr{i}(1,3),sprintf('%d',i));
        end
        
    end
end