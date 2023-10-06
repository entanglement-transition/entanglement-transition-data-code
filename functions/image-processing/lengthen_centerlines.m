function new_cl_list = lengthen_centerlines(zstack,cl_list)

num_rods = numel(cl_list);
N = 1000;

new_cl_list = cell(num_rods,1);
for i = 1:num_rods
    cl = cl_list{i};
    fr = fit_rod(cl');

    delta_th = (max(fr.philist) - min(fr.philist))/(N/4);
    min_th = min(fr.philist) - delta_th*(N/4);
    max_th = max(fr.philist) + delta_th*(N/4);
    ex_philist = linspace(min_th,max_th,N)';

    u1 = fr.u(:,1);
    u2 = fr.u(:,2);
    xpt2 = fr.xm+fr.r*cos(ex_philist)';
    ypt2 = fr.ym+fr.r*sin(ex_philist)';

    extended_rod = (fr.cen + xpt2.*u1 + ypt2.*u2)';
    extended_rod = round(unique(extended_rod,'rows','stable'));
    extended_rod = cut_to_matrix_size(size(zstack),extended_rod);

    ind = sub2ind2(size(zstack),extended_rod);
    true_ind = zstack(ind);

    new_cl_list{i} = extended_rod(true_ind,:);
end

end