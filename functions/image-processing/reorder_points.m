function [xyz2, I] = reorder_points(xyz,th,start_index)

% what about initial_point?
N = size(xyz,1);

I(1) = start_index;
% xyz = setdiff(xyz,xyz(1,:),'rows');

for i = 2:N

    dist = sqrt(sum( (xyz - xyz(I(end),:)).^2 , 2));        
    [d_sorted,I_sorted] = sort(dist);        
    
    j = 2;
    while 1
        
        if ~ismember(I_sorted(j),I)
            
            if dist(I_sorted(j)) > th
                break
            end
                
            I(end+1) = I_sorted(j);
            break            
        else
            j = j + 1;
        end
    end
    
end
    


xyz2 = xyz(I,:);

% xy = d;
% sortedxy = sortrows(xy,2); % sort Y-values relative to X-values in ascending order
% [~,idu] = unique(sortedxy(:,2)); % get the index of unique Y-values
% for ii = 2:2:length(idu)-1
% 
%     sortedxy(idu(ii-1):idu(ii)-1,1) = sort(sortedxy(idu(ii-1):idu(ii)-1,1));  % sort X-values 
%     % relative to odd indexed unique Y values in ascending order
%     sortedxy(idu(ii):idu(ii+1)-1,1) = sort(sortedxy(idu(ii):idu(ii+1)-1,1),'descend'); % sort X-values
%     % relative to even indexed unique Y values in descending order
% end

end