classdef FieldArray
    
    properties
        fields
        
        cx
        cy
        cz
        R
        h
                
        g
    end
    
    methods        
        function obj = FieldArray(field_struct,field_list)
            % TO DO: 3d field as FieldAnalysis? or more primitive ones...
            % like?
            if ismember('c',field_list)
                obj.fields = containers.Map({'e','phi','n','s','c'},...
                    {field_struct.e,field_struct.phi,field_struct.n,field_struct.s,field_struct.c});
            else
                obj.fields = containers.Map({'e','phi','n','s'},...
                    {field_struct.e,field_struct.r,field_struct.n,field_struct.s});
            end
            
            obj.cx = field_struct.cx;
            obj.cy = field_struct.cy;
            obj.cz = field_struct.cz;
            obj.R = field_struct.R;
            obj.h = field_struct.h;
            
        end
        
        function [cropped_fields,rmask,crop] = crop(obj,varargin)
            if nargin == 2
                cutoff = varargin{1};
                zcutoff = 0;
            elseif nargin == 3
                cutoff = varargin{1};
                zcutoff = varargin{2};
            elseif nargin == 1
                cutoff = 0.9;
                zcutoff = 0;
            end
            
            % based on density field
            rho = obj.fields('r');
            horizontal_projection = rot90(squeeze(max(rho,[],1)));
            vertical_projection = max(rho,[],3);
            tmp = squeeze( sum(sum(rho,1),2) );
            tmp2 = tmp> mean(tmp(:));
            
            zcut = max(find(tmp2))*0.9; %0.8 for 0.76 mm
            zcut_low = min(find(tmp2))*zcutoff; % 1.2 for 0.76 mm
            zcut_low = max([zcut_low 1]);
            
            bw = vertical_projection>0;
            bw2 = padarray(bw,[10 10]);
            
            %
            %             figure
            %             imagesc(bw2);
            
            [a,b,c] = size(rho);
            [centers, radii, metric] = imfindcircles(bw2,round([b/3 2*b/3]),'Sensitivity',0.85);
            cx = centers(2)-10;
            cy = centers(1)-10;
            
            R = radii*cutoff; % 0.85 for 0.76 mm
            %     zcut = round(size(rho,3)*90/113);
            rcut = R;
            
            rmask = get_ring_mask(size(rho),0,rcut);
            crop = [1 1 zcut_low size(rho,1)-1 size(rho,1)-1 zcut-zcut_low-1];
            % visual check?
            
            cropped_fields = containers.Map(keys(obj.fields),values(obj.fields));
            for k = keys(obj.fields)
                fd = cropped_fields(k{1});
                fd(~rmask) = NaN;
                fd = imcrop3(fd,crop);
                cropped_fields(k{1}) = fd;
            end
            
            fig_on = 0;
            if fig_on
                
                pos_hori = [cy-rcut c-zcut 2*rcut zcut-zcut_low];
                pos_vert = [cy-rcut cx-rcut 2*rcut 2*rcut];
                
                f = figure('position',[1885 438 1027 980]);
                pos1 = [0.1 0.4 0.6 0.2];
                subplot('position',pos1)
                imagesc(horizontal_projection)
                rectangle('Position',pos_hori,'edgecolor','r','linewidth',1)
                rectangle('Position',round([b/2 c/2 obj.w(2)/(obj.w(2)-obj.o(2)) obj.w(3)/(obj.w(3)-obj.o(3))]),'edgecolor','r','linewidth',1)
                axis equal
                axis off
                xlim_tmp = get(gca,'XLim');
                ylim_tmp = get(gca,'YLim');
                
                pos2 = [0.1 0.1 0.6 0.22];
                subplot('position',pos2)
                plot(mean(horizontal_projection,1))
                xlim(xlim_tmp)
                ylim_tmp2 = get(gca,'YLim');
                hold on
                plot([pos_vert(1),pos_vert(1)],ylim_tmp2,'r-')
                plot([pos_vert(1)+pos_vert(3),pos_vert(1)+pos_vert(3)],ylim_tmp2,'r-')
                xlabel('$x$ (px)')
                ylabel('Max. proj.')
                
                pos3 = [0.7 0.4 0.2 0.2];
                subplot('position',pos3)
                tmp = mean(horizontal_projection,2);
                plot(flip(tmp),1:size(obj.fields('r'),3))
                ylim(ylim_tmp)
                xlim_tmp = get(gca,'XLim');
                hold on
                plot(xlim_tmp,[zcut, zcut],'r-')
                xlabel('Max. proj.')
                ylabel('$z$ (px)')
                
                pos4 = [0.1 0.65 0.6 0.344];
                subplot('position',pos4)
                imagesc(vertical_projection)
                rectangle('Position',pos_vert,'Curvature',[1 1],'edgecolor','r','linewidth',1)
                rectangle('Position',round([a/2 b/2 obj.w(1)/(obj.w(1)-obj.o(1)) obj.w(2)/(obj.w(2)-obj.o(2))]),'edgecolor','r','linewidth',1)
                axis equal
                axis off
            end
        end
        
        function avg_fields = average_fields(obj_array)
            N = numel(obj_array);
            clear sz
            for i_scan = 1:N
                davg_fd = mean(obj_array(i_scan).fields('e') ,3,'omitnan');
                havg_fd = rot90(squeeze(mean(obj_array(i_scan).fields('e'),2,'omitnan')));
                
                sz(i_scan).dsz = size(davg_fd);
                sz(i_scan).hsz = size(havg_fd);
            end            
            hsz_list = vertcat(sz.hsz);
            HSZ = max(hsz_list(:,1));
            for i_scan = 1:N                
                da = containers.Map(keys(obj_array(i_scan).fields),values(obj_array(i_scan).fields));
                ha = containers.Map(keys(obj_array(i_scan).fields),values(obj_array(i_scan).fields));

                for k = keys(da)
                    davg_fd = mean(obj_array(i_scan).fields(k{1}) ,3,'omitnan');
                    havg_fd = rot90(squeeze(mean(obj_array(i_scan).fields(k{1}),2,'omitnan')));
                    
                    num_added = HSZ-size(havg_fd,1);
                    tmp = padarray(havg_fd,num_added,'pre');

                    da(k{1}) = davg_fd;
                    ha(k{1}) = tmp;
%                     size(obj_array(i_scan).fields(k{1}))
                end                
                avg_fields(i_scan).da = da;
                avg_fields(i_scan).ha = ha;                
            end

        end

        function plot_horizontal_section(obj,s)
            
        end
        function calculate_variance(obj)
            
        end
        function [out,out2] = get_representative_values(obj_array,s,p)
            % h: handle
            N = numel(obj_array);
            out = zeros(1,N);
            out2 = zeros(1,N);
            for i = 1:N
                obj = obj_array(i);
                target = obj.fields(s);
                out(i) = yj.p_norm(target,p);
                out2(i) = std(target(:))/numel(target(:));
%                 nf = obj.fields('n');
%                 out(i) = mean(target(nf>1),'all');
            end
        end
        
        function out = plot_statistics(obj_array,s)
            N = numel(obj_array);
            out = zeros(1,N);
            for i = 1:N
                obj = obj_array(i);
                target = obj.fields(s);
                histogram(target);
                hold on
            end
        end
        
%         function out = (obj_array)
%             % h: handle
%             N = numel(obj_array);
%             out = zeros(1,N);
%             for i = 1:N
%                 
%             end
%         end
        
        
    end
    
    
end


function [crop, rmask] = get_ROI(map,varargin)

if nargin == 2
    cutoff = varargin{1};
    zcutoff = 0;
elseif nargin == 3
    cutoff = varargin{1};
    zcutoff = varargin{2};
elseif nargin == 1
    cutoff = 0.9;
    zcutoff = 0;
end

horizontal_projection = rot90(squeeze(max(map.r,[],1)));
vertical_projection = max(map.r,[],3);

%     figure
%     tmp = squeeze( sum(sum(map.r,1),2) );
%     plot(flip(tmp),1:size(map.r,3),'o-')
%
%     hold on
%     tmp2 = tmp>0;
%     plot(tmp2,1:size(map.r,3),'o-')

tmp = squeeze( sum(sum(map.r,1),2) );
tmp2 = tmp> mean(tmp(:));

zcut = max(find(tmp2))*0.9; %0.8 for 0.76 mm
zcut_low = min(find(tmp2))*zcutoff; % 1.2 for 0.76 mm
zcut_low = max([zcut_low 1]);

bw = vertical_projection>0;
bw2 = padarray(bw,[10 10]);

%
%     figure
%     imagesc(bw2);

[a,b,c] = size(map.r);
[centers, radii, metric] = imfindcircles(bw2,round([b/3 2*b/3]),'Sensitivity',0.85);
cx = centers(2)-10;
cy = centers(1)-10;

R = radii*cutoff; % 0.85 for 0.76 mm
%     zcut = round(size(map.r,3)*90/113);
rcut = R;

rmask = get_ring_mask(size(map.r),0,rcut);
crop = [1 1 zcut_low size(map.r,1)-1 size(map.r,1)-1 zcut-zcut_low-1];

cropped = map.r;
cropped(~rmask) = NaN;
cropped = imcrop3(cropped,crop);

pos_hori = [cy-rcut c-zcut 2*rcut zcut-zcut_low];
pos_vert = [cy-rcut cx-rcut 2*rcut 2*rcut];

f = figure('position',[1885 438 1027 980]);
pos1 = [0.1 0.4 0.6 0.2];
subplot('position',pos1)
imagesc(horizontal_projection)
rectangle('Position',pos_hori,'edgecolor','r','linewidth',1)
rectangle('Position',round([b/2 c/2 obj.w(2)/(map.w(2)-map.o(2)) map.w(3)/(map.w(3)-map.o(3))]),'edgecolor','r','linewidth',1)
axis equal
axis off
xlim_tmp = get(gca,'XLim');
ylim_tmp = get(gca,'YLim');

pos2 = [0.1 0.1 0.6 0.22];
subplot('position',pos2)
plot(mean(horizontal_projection,1))
xlim(xlim_tmp)
ylim_tmp2 = get(gca,'YLim');
hold on
plot([pos_vert(1),pos_vert(1)],ylim_tmp2,'r-')
plot([pos_vert(1)+pos_vert(3),pos_vert(1)+pos_vert(3)],ylim_tmp2,'r-')
xlabel('$x$ (px)')
ylabel('Max. proj.')

pos3 = [0.7 0.4 0.2 0.2];
subplot('position',pos3)
tmp = mean(horizontal_projection,2);
plot(flip(tmp),1:size(map.r,3))
ylim(ylim_tmp)
xlim_tmp = get(gca,'XLim');
hold on
plot(xlim_tmp,[zcut, zcut],'r-')
xlabel('Max. proj.')
ylabel('$z$ (px)')

pos4 = [0.1 0.65 0.6 0.344];
subplot('position',pos4)
imagesc(vertical_projection)
rectangle('Position',pos_vert,'Curvature',[1 1],'edgecolor','r','linewidth',1)
rectangle('Position',round([a/2 b/2 map.w(1)/(map.w(1)-map.o(1)) map.w(2)/(map.w(2)-map.o(2))]),'edgecolor','r','linewidth',1)
axis equal
axis off


end

function fd = get_fd(fields,field_type)

switch field_type
    case 'e'
        fd = fields.e;
    case 'r'
        fd = fields.r;
    case 'n'
        fd = fields.n;
    case 's'
        fd = fields.s;
    case 'ne1' % just in case
        fd = fields.ne1;
    case 'ne2'
        fd = fields.ne2;
end

end

function [crop, rmask] = get_ROI_for_7(map)

horizontal_projection = rot90(squeeze(max(map.r,[],1)));
vertical_projection = max(map.r,[],3);

%     figure
%     tmp = squeeze( sum(sum(map.r,1),2) );
%     plot(flip(tmp),1:size(map.r,3),'o-')
%
%     hold on
%     tmp2 = tmp>0;
%     plot(tmp2,1:size(map.r,3),'o-')

tmp = squeeze( sum(sum(map.r,1),2) );
tmp2 = tmp> mean(tmp(:));
zcut = max(find(tmp2))*0.9; %0.8 for 0.76 mm
zcut_low = min(find(tmp2))*1.1; % 1.2 for 0.76 mm
zcut_low = max([zcut_low 1]);

bw = vertical_projection>0;
bw2 = padarray(bw,[10 10]);

%
%     figure
%     imagesc(bw2);

[a,b,c] = size(map.r);
[centers, radii, metric] = imfindcircles(bw2,round([b/3 2*b/3]),'Sensitivity',0.9);
cx = centers(2)-10;
cy = centers(1)-10;

R = radii*0.9; % 0.85 for 0.76 mm
%     zcut = round(size(map.r,3)*90/113);
rcut = R;

rmask = get_ring_mask(size(map.r),0,rcut);
crop = [1 1 zcut_low size(map.r,1)-1 size(map.r,1)-1 zcut-zcut_low-1];

cropped = map.r;
cropped(~rmask) = NaN;
cropped = imcrop3(cropped,crop);

pos_hori = [cy-rcut c-zcut 2*rcut zcut-zcut_low];
pos_vert = [cy-rcut cx-rcut 2*rcut 2*rcut];

f = figure('position',[1885 438 1027 980]);
pos1 = [0.1 0.4 0.6 0.2];
subplot('position',pos1)
imagesc(horizontal_projection)
rectangle('Position',pos_hori,'edgecolor','r','linewidth',1)
rectangle('Position',round([b/2 c/2 map.w(2)/(map.w(2)-map.o(2)) map.w(3)/(map.w(3)-map.o(3))]),'edgecolor','r','linewidth',1)
axis equal
axis off
xlim_tmp = get(gca,'XLim');
ylim_tmp = get(gca,'YLim');

pos2 = [0.1 0.1 0.6 0.22];
subplot('position',pos2)
plot(mean(horizontal_projection,1))
xlim(xlim_tmp)
ylim_tmp2 = get(gca,'YLim');
hold on
plot([pos_vert(1),pos_vert(1)],ylim_tmp2,'r-')
plot([pos_vert(1)+pos_vert(3),pos_vert(1)+pos_vert(3)],ylim_tmp2,'r-')
xlabel('$x$ (px)')
ylabel('Max. proj.')

pos3 = [0.7 0.4 0.2 0.2];
subplot('position',pos3)
tmp = mean(horizontal_projection,2);
plot(flip(tmp),1:size(map.r,3))
ylim(ylim_tmp)
xlim_tmp = get(gca,'XLim');
hold on
plot(xlim_tmp,[zcut, zcut],'r-')
xlabel('Max. proj.')
ylabel('$z$ (px)')

pos4 = [0.1 0.65 0.6 0.344];
subplot('position',pos4)
imagesc(vertical_projection)
rectangle('Position',pos_vert,'Curvature',[1 1],'edgecolor','r','linewidth',1)
rectangle('Position',round([a/2 b/2 map.w(1)/(map.w(1)-map.o(1)) map.w(2)/(map.w(2)-map.o(2))]),'edgecolor','r','linewidth',1)
axis equal
axis off


end

function mask_3d = get_ring_mask(sz,ri,ro,varargin)

if nargin == 4
    dv =varargin{1};
elseif nargin == 3
    dv = [0 0 0];
end

a = sz(1);
b = sz(2);
c = sz(3);

cx = round(a/2) + dv(1);
cy = round(b/2) + dv(2);
cz = round(c/2) + dv(3);

[X,Y] = meshgrid( (1:a) - cx, (1:b) - cy);
R = sqrt(X.^2 + Y.^2);

mask = (ri - 1 < R & R < ro + 1);
mask_3d = repmat(mask,[1 1 c]);

end
