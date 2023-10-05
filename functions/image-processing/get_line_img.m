function limg = get_line_img(image_size,points)    
    assert(size(points,2) == 3);
    limg = zeros(image_size,'logical');    
    I = sub2ind2(image_size,points);
    limg(I) = 1;
end