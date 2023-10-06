function zstack2 = clear_zstack(zstack2,volume_threshold)

cc = bwconncomp(zstack2);
stat = regionprops3(cc,'Volume');
I_small = stat.Volume < volume_threshold;
new_pixel_idx_list = cc.PixelIdxList;
new_pixel_idx_list(I_small) = [];
num_obj = numel(new_pixel_idx_list);
cc2 = cc;
cc2.NumObjects = num_obj;
cc2.PixelIdxList = new_pixel_idx_list;
zstack2 = binarymatrix(cc2);

end