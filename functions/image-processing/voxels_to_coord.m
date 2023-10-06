function [cl_list,error_list] = voxels_to_coord(zstack,N)

cc = bwconncomp(zstack);

ind_list = cc.PixelIdxList;
num_elements = cellfun(@numel, ind_list);
ind_list(num_elements<N) = [];


num_cc = numel(ind_list);
error_list = zeros(num_cc,1);
cl_list = cell(num_cc,1);

tstart = tic;
for i = 1:num_cc
    ind = ind_list{i};
    cl = ind2sub2(size(zstack),ind);
    
    [cen,ori,slist] = get_line_coord(cl);
    [~,I] = sort(slist);
    cl = cl(I,:);
    cl_list{i} = cl;
    
    dlist = rwnorm(cl - (cen + ori.*slist));
    error_list(i) = mean(dlist);
     if mod(i,100) == 0
        fprintf('%d / %d processed. %.2f sec elapsed.\n',i,num_cc,toc(tstart));
    end
end



end