function points = cut_to_matrix_size(matrix_size,points)
    assert(size(points,2) == 3)
    points(find(any(points <= 0,2)),:) = [];
    points(find(any(points >= matrix_size,2)),:) = [];
    
end