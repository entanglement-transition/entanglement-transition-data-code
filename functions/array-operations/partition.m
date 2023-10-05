function parts = partition(collection,k,enemies)
    if length(collection) == 1
        parts = {num2cell(collection)};
        return;
    end
    
    first = collection(1);
    rest = collection(2:end);
    
    % Recursive call
    rest_parts = partition(rest,k,enemies);

    parts = {};

    % For each smaller partition...
    for i = 1:length(rest_parts)
        rest = rest_parts{i};

        % Insert `first` in each of the subpartition's subsets
        for j = 1:length(rest)
            new_subset = [first, rest{j}];
            new_part = rest;
            new_part{j} = new_subset;
            
            if length(new_part) <= k % Check the number of subsets
                parts{end+1} = new_part;
            end
%             parts{end+1} = new_part;
        end

        % Put `first` in its own subset 
        new_part = [[first], rest];
        if length(new_part) <= k % Check the number of subsets
            parts{end+1} = new_part;
        end
%         parts{end+1} = [[first], rest];
    end

    to_delete = zeros(numel(parts),1,'logical');
    for i = 1:numel(parts)
        part_i = parts{i};
        I = cellfun(@(x) any(all(ismember(enemies,x),2),'all'),part_i);
        
        if any(I)
            to_delete(i) = 1;
        end        
    end
    parts(to_delete) = [];
    
end
