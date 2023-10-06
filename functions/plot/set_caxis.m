function set_caxis(h,map,cmp,VAxis,varargin)
    
    if nargin == 4 | (nargin == 4 & strcmp(varargin{1},'lin'))
        % linear by default
        newscale = linspace(min(map(:)) - min(VAxis(:)),...
        max(map(:)) - min(VAxis(:)), size(cmp, 1))/diff(VAxis);
        newscale(newscale < 0) = 0;
        newscale(newscale > 1) = 1;
        cmp2 = interp1(linspace(0, 1, size(cmp, 1)), cmp, newscale);

        h.Colormap = cmp2;
    elseif strcmp(varargin{1},'log')
        newscale = logspace( log10(min(map(:)) - min(VAxis(:))),...
        log10(max(map(:)) - min(VAxis(:))), size(cmp, 1))/diff(VAxis);
    
        newscale(newscale < 0) = 0;
        newscale(newscale > 1) = 1;
        cmp2 = interp1(linspace(0, 1, size(cmp, 1)), cmp, newscale);    

        h.Colormap = cmp2;
    end
end