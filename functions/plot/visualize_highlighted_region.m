function [h,lifted,flash] = visualize_highlighted_region(stack,highlight)
    highlight = unique(round(highlight),'rows','stable');
    highlight = cut_to_matrix_size(size(stack),highlight);
    cuboid = get_bbox(size(stack),highlight);

    lifted = imcrop3(stack,cuboid);
    flash = get_line_img(size(lifted),highlight - cuboid([2 1 3]) + [1 1 1]);


    flash = imdilate(flash,strel('sphere',2));
    A = lifted + 2*flash;


    num_labels = numel(unique(A(:)));

    label_opacity = ones(num_labels,1)*0.04;
    label_opacity(1) = 0;
    label_opacity(end) = 0.8;
    label_color = ones(num_labels,3)*0.1;
    label_color(end-1,:) = [0 0 1];
    label_color(3,:) = [0 1 0];
    label_color(end,:) = [1 1 0];

    
    h = labelvolshow(A);
    h.LabelOpacity = label_opacity;
    h.LabelColor = label_color;
    h.BackgroundColor = 'w';

end