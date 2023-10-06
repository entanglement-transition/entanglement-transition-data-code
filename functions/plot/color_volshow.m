function h = color_volshow(vol,cmap,caxis)
    h = volshow(vol,'renderer','maximumintensityprojection','backgroundcolor',[0 0 0],'colormap',jet);    
    if nargin > 2
        set_caxis(h,vol,cmap,caxis);
    end
end
