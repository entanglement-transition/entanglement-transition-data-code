function plot_projected_fields(field)
% plot averaged fields in horizontally and vertically
% not projections in a rigorous terms, but...

h = set_figure(15,15);
tiledlayout(num_samples,4);
for i_sample = 1:num_samples
    sample_id = string(samples(i_sample));
    FA = FA_all(sample_id);
    
    zall = zeros(1,4);
    for i = 1:4
        [z1,~] = size(FA(i).average_fields.ha(fstr));
        zall(i) = z1;
    end
    zmax = max(zall);
    
    
    
    for i_strain = 1:4
        ef = FA(i_strain).average_fields.ha(fstr);
        
        [z1,~] = size(FA(i_strain).average_fields.ha(fstr));
        deficit = zmax - z1;
        ef = padarray(ef,[deficit,1],'pre');
        
        if strcmp(fstr,'e')
            ef = ef/2;
        end
        
        if strcmp(fstr,'c')
            rf = FA(i_strain).average_fields.ha('r');
            load(path_all(sample_id).contact_path{i_scan});
            num_all_contacts = size(contact_info_all,1);
            total_vol = prod(stack_size_list{i_sample}(i_scan,:));
            norm_fac = 4/3*pi*(FA(i_strain).R)^3/total_vol*num_all_contacts;
            ef = ef/norm_fac;
        end
        
        nexttile
        imagesc(ef);hold on;
        colormap(coolwarm)
        caxis([0,max_field(i_fstr)]);
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])
        %         axis off
        
        if i_strain == 1
            ylabel(sprintf('$\\alpha = %d$',alpha_list(i_sample)),'interpreter','latex');
        end
        
        if i_sample == 1
            title(sprintf('$\\epsilon = %.2f$',eps_list(i_strain)),'interpreter','latex');
        end
    end
    
end
colorbar
%
image_name = sprintf('all_fields_ha_%s.png',fstr);
% panel_path = fullfile(dropbox_location,'ResearchFigures/[Entanglement]Supplementary/panels');

end