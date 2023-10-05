function d = get_root_dir()
% Get root directory of toolbox
    sp = split(path,":");
    for k=1:length(sp)
        p = cell2mat(sp(k));
        if(p(end-6:end) == "pov_lab")
            d = p;
            return;
        end
    end
end