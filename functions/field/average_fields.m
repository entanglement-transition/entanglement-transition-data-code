function [vertically_averaged_field,horizontally_averaged_field] = average_fields(all_fields)
% all_fields: a structure variable



for i_field = 1:5
    
    % to do: adapt to the number of fields in the input structure?
    field = choose_field(all_fields,i_field);
    
    davg_fd = mean(obj_array(i_scan).fields(k{1}) ,3,'omitnan');
    havg_fd = rot90(squeeze(mean(obj_array(i_scan).fields(k{1}),2,'omitnan')));
    
    num_added = HSZ-size(havg_fd,1);
    tmp = padarray(havg_fd,num_added,'pre');
    
    da(k{1}) = davg_fd;
    ha(k{1}) = tmp;
    
end



end