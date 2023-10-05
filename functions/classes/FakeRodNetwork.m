classdef FakeRodNetwork < handle
    properties
        main_rod % CircularRod instance
        contacting_rods % Array of CircularRod instances
        
        all_rods
        contact_table
    end
    
    methods
        function obj = FakeRodNetwork(main_rod)
            if nargin > 0
                obj.main_rod = main_rod;                
                obj.all_rods = [obj.main_rod];
            end
        end
        
        function found1_or_not0 = create_contacting_rods_randomly(obj,num_contacts)
            R0 = obj.main_rod.rod_radius;
            L0 = obj.main_rod.rod_length;
            phi_min0 = obj.main_rod.phi_min;
            phi_max0 = obj.main_rod.phi_max;
            
            % random sampling
            max_outeriter = 200;
            max_inneriter = 50;
            outeriter = 1;
            while outeriter < max_outeriter
                phi_contacts = (phi_max0 - phi_min0)*rand(num_contacts,1) + phi_min0;
                theta_contacts = 2*pi*rand(num_contacts,1);
                contact_foots = obj.main_rod.query_points(phi_contacts);
                directors = obj.main_rod.query_directors(phi_contacts);
                
                contact_points = contact_foots ...
                    + R0*cos(theta_contacts).*directors(:,1:3) ...
                    + R0*sin(theta_contacts).*directors(:,4:6);
                contact_foots_j = zeros(size(contact_foots));
                contacting_rods(1:num_contacts) = CircularRod;
                Q_i = zeros(num_contacts,9);
                Q_j = zeros(num_contacts,9);
                for i = 1:num_contacts
                    v_ij = contact_points(i,:) - contact_foots(i,:);
                    contact_foots_j(i,:) = contact_foots(i,:) + 2*v_ij; % contact foot at rod j
                    
                    inneriter = 1;
                    while inneriter < max_inneriter
                        rc = uniform_log_random(1,3,4);
                        while 1
                            d3_at_foot = cross(v_ij,sphere_pick(1)); % tangent at the contact point
                            d3_at_foot = d3_at_foot/norm(d3_at_foot);
                            if abs( dot(d3_at_foot, v_ij/norm(v_ij)) ) < 1e-15
                                break;
                            end
                        end
                        
                        random_sign = sign(rand(1) - 0.5);
                        d1_at_foot = random_sign*v_ij/norm(v_ij);
                        %                         d1_at_foot = -v_ij/norm(v_ij);
                        d2_at_foot = cross(d3_at_foot,d1_at_foot);
                        d2_at_foot = d2_at_foot/norm(d2_at_foot);
                        triad = [d1_at_foot,d2_at_foot,d3_at_foot];
                        phi_max = rand(1)*L0/rc;
                        phi_min = phi_max - L0/rc;
                        
                        n_ij = -v_ij/norm(v_ij);
                        b_ij = cross(directors(i,7:9),n_ij);
                        t_ij = cross(n_ij,b_ij);
                        
                        n_ji = -n_ij;
                        b_ji = cross(d3_at_foot,n_ji);
                        t_ji = cross(n_ji,b_ji);
                        
                        Q_i(i,:) = [n_ij,b_ij,t_ij];
                        Q_j(i,:) = [n_ji,b_ji,t_ji];
                        
                        candidate = CircularRod(contact_foots_j(i,:),triad,rc,R0,L0,phi_min,phi_max);
                        if get_distance(candidate,obj.main_rod) < 2*R0;
                            ; % go on
                        elseif is_overlap(candidate,contacting_rods)
                            ; % go on
                        else
                            contacting_rods(i) = candidate;
                            break;
                        end
                        
                        inneriter = inneriter + 1;
                        
                        if 0
                            close all;
                            pts = obj.main_rod.query_points(phi_contacts);
                            plot(candidate);hold on;
                            plot(obj.main_rod);
                            plot(pts,'ro');
                            ;
                        end
                    end
                    
                    if inneriter == max_inneriter
                        break;
                    end
                end
                
                if 0
                    %                 pts = obj.main_rod.get_discrete_points(100);
                    
                    close all;
                    plot(obj.main_rod);hold on;
                    plot(candidate);hold on;
                    plot3v(contact_foots_j,'ko');
                    
                    plot3v(contact_foots,'ro');
                    plot3v(contact_points,'bo');
                    
                    ;axis equal;
                end
                
                outeriter = outeriter + 1;
                ind_nonempty = ~isempty(contacting_rods);
                if all(ind_nonempty)
                    break;
                end
                disp(sprintf('Finding contacting rods...\t Iteration: %d',outeriter));
            end
            
            found1_or_not0 = all(~isempty(contacting_rods));
            obj.contacting_rods = contacting_rods;
            obj.all_rods = [obj.main_rod,obj.contacting_rods];
            
            i = ones(num_contacts,1);
            j = (2:num_contacts+1)';
            contact_table = table(i,j,contact_points,contact_foots,contact_foots_j,...
                Q_i,Q_j,phi_contacts,theta_contacts);
            contact_table = renamevars(contact_table,...
                contact_table.Properties.VariableNames, ...
                ["i","j","p_i","s_i","s_j","Q_i","Q_j","phi_c","theta_c"]);
            obj.contact_table = contact_table;
            
        end
        
        function [new_FRN,found1_or_not0] = perturb_randomly(obj)
            % m-rod
            A_translation = 3; % rigid body translation, in pixel, normal dist
            A_rolling = 1; % contact rolling, in degree
            
            % both m-rod and c-rods
            A_curvature = 0.05; % ratio of increased curvature, uniform dist
            A_rotation = 1; % variance of rigid body rotation, in degree, normal dist
            A_sliding = 3; % contact sliding, in pixel
            
            % c-rods
            A_rotation_c = 1;
            
            % perturb main rod first
            R0 = obj.main_rod.rod_radius;
            L0 = obj.main_rod.rod_length;
            
            rc0 = obj.main_rod.radius_curvature;
            triad0 = obj.main_rod.plane_triad;
            origin0 = obj.main_rod.plane_origin;
            
            new_rc = rc0 * (1 + rand(1)*A_curvature);
            rotation_axis = sphere_pick(1); % uniform over sphere
            K = [0 -rotation_axis(3) rotation_axis(2) ; rotation_axis(3) 0 -rotation_axis(1) ; -rotation_axis(2) rotation_axis(1) 0 ];
            th = A_rotation*randn(1)*pi/180; % rotation mag.: normal distribution
            Q = eye(3,3) + sin(th)*K + (1-cos(th))*K^2;
            
            new_triad = Q*reshape(triad0,[3,3]);
            new_triad = reshape(new_triad,[1,9]);
            translation = A_translation*randn(1)*sphere_pick(1); % normal dist
            new_origin = obj.main_rod.plane_origin + translation;
            phi_min0 = obj.main_rod.phi_min;
            phi_max0 = obj.main_rod.phi_max;
            
            new_phi_min = -L0/new_rc/2 + (phi_min0 + phi_max0)/2; % should preserve arclength
            new_phi_max = L0/new_rc/2 + (phi_min0 + phi_max0)/2;
            
            new_main_rod = CircularRod(new_origin,new_triad,new_rc,R0,L0,new_phi_min,new_phi_max);
            new_FRN = FakeRodNetwork(new_main_rod);
            
            % now perturb contacting rods - tricky part
            % new phi and theta should be around the former ones
            phi_c = obj.contact_table.phi_c;
            theta_c = obj.contact_table.theta_c;
            num_contacts = numel(phi_c);
            
            max_outeriter = 200;
            max_inneriter = 50;
            outeriter = 1;
            while outeriter < max_outeriter
                % tricky part
                inneriter = 1;
                while inneriter < max_inneriter
                    d = A_sliding*randn(num_contacts,1)*(phi_max0 - phi_min0)/L0; % normal dist
                    phi_contacts = phi_c * (new_phi_max-new_phi_min)/(phi_max0-phi_min0) + d;
                    
                    I_valid = ((phi_contacts <= new_phi_max) & ...
                        (phi_contacts >= new_phi_min));
                    ;
                    if all( I_valid )
                        break;
                    end
                    phi_contacts(~I_valid) = phi_c(~I_valid);
                    inneriter = inneriter + 1;
                end
                theta_contacts = A_rolling*randn(1)*pi/180 + theta_c;
                % perturbed
                contact_foots = new_main_rod.query_points(phi_contacts);
                directors = new_main_rod.query_directors(phi_contacts);
                
                %                 rot = [cos(theta_contacts),sin(theta_contacts);...
                %                     -sin(theta_contacts),cos(theta_contacts)];
                
                contact_points = contact_foots ...
                    + R0*cos(theta_contacts).*directors(:,1:3) ...
                    + R0*sin(theta_contacts).*directors(:,4:6);
                
                %                 dot(contact_points(1,:) - contact_foots(1,:), directors(7:9))
                %                 dot(directors(1,1:3), directors(1,7:9))
                %                 sum(directors(:,1:3).*directors(:,7:9),2)
                %                 sum(( contact_points - contact_foots ).*directors(:,7:9),2)
                %                 sum(( contact_points - contact_foots ).*directors(:,7:9),2)
                %                 directors(:,1:3)
                %                 directors(:,4:6)
                
                contact_foots_j = zeros(size(contact_foots));
                contacting_rods(1:num_contacts) = CircularRod;
                Q_i = zeros(num_contacts,9);
                Q_j = zeros(num_contacts,9);
                for i = 1:num_contacts
                    v_ij = contact_points(i,:) - contact_foots(i,:);
                    contact_foots_j(i,:) = contact_foots(i,:) + 2*v_ij; % contact foot at rod j
                    
                    inneriter = 1;
                    while inneriter < max_inneriter
                        rc = obj.contacting_rods(i).radius_curvature * (1 + rand(1)*A_curvature);
                        
                        rotation_axis = sphere_pick(1); % uniform over sphere
                        K = [0 -rotation_axis(3) rotation_axis(2) ; rotation_axis(3) 0 -rotation_axis(1) ; -rotation_axis(2) rotation_axis(1) 0 ];
                        th = A_rotation_c*randn(1)*pi/180; % rotation mag.: normal distribution
                        Q = eye(3,3) + sin(th)*K + (1-cos(th))*K^2;
                        old_triad = obj.contacting_rods(i).plane_triad;
                        
                        contact_foots0 = obj.main_rod.query_points(phi_c(i));
                        directors0 = obj.main_rod.query_directors(phi_c(i));
                        contact_points0 = contact_foots0 ...
                            + R0*cos(theta_c(i)).*directors0(:,1:3) ...
                            + R0*sin(theta_c(i)).*directors0(:,4:6);
                        v_ij0 = contact_points0 - obj.main_rod.query_points(phi_c(i));
                        
                        Q_m1 = reshape(obj.main_rod.query_directors(phi_c(i)),[3,3])';
                        Q_m2 = reshape(directors(i,:),[3,3])';
                        Q_c1 = reshape(obj.contacting_rods(i).query_directors(0),[3,3])';
                        Q_c2 = Q_m2*Q_m1'*Q_c1; % natural triad - not necessarily?
                        
                        d3_temp = Q_c2(3,:);
                        
                        rotation_axis = sphere_pick(1); % uniform over sphere
                        K = [0 -rotation_axis(3) rotation_axis(2) ; rotation_axis(3) 0 -rotation_axis(1) ; -rotation_axis(2) rotation_axis(1) 0 ];
                        th = A_rotation*randn(1)*pi/180; % rotation mag.: normal distribution
                        Q = eye(3,3) + sin(th)*K + (1-cos(th))*K^2;
                        d3_temp = (Q*d3_temp')';
                        
                        d1_at_foot = sign(dot(v_ij0,old_triad(1:3)))*v_ij/norm(v_ij);
                        d2_at_foot = cross(d3_temp,d1_at_foot);
                        d2_at_foot = d2_at_foot/norm(d2_at_foot);
                        d3_at_foot = cross(d1_at_foot,d2_at_foot);
                        d3_at_foot = d3_at_foot/norm(d3_at_foot);
                        triad = [d1_at_foot,d2_at_foot,d3_at_foot];
                        
                        phi_min = -L0/rc/2 + (obj.contacting_rods(i).phi_min + obj.contacting_rods(i).phi_max)/2; % should preserve arclength
                        phi_max = L0/rc/2 + (obj.contacting_rods(i).phi_min + obj.contacting_rods(i).phi_max)/2;
                        
                        while 1
                            d = A_sliding*randn(1)*(obj.contacting_rods(i).phi_max -...
                                obj.contacting_rods(i).phi_min)/L0; % normal dist
                            phi_min_p = phi_min + d;
                            phi_max_p = phi_max + d;
                            
                            if (phi_min_p <= 0) & (phi_max_p >= 0)
                                break;
                            end
                        end
                        
                        %                         (phi_max_p - phi_min_p)*rc
                        n_ij = -v_ij/norm(v_ij);
                        b_ij = cross(directors(i,7:9),n_ij);
                        t_ij = cross(n_ij,b_ij);
                        
                        n_ji = -n_ij;
                        b_ji = cross(d3_at_foot,n_ji);
                        t_ji = cross(n_ji,b_ji);
                        
                        Q_i(i,:) = [n_ij,b_ij,t_ij];
                        Q_j(i,:) = [n_ji,b_ji,t_ji];
                        candidate = CircularRod(contact_foots_j(i,:),triad,rc,R0,L0,phi_min_p,phi_max_p);
                        %                         while abs(get_distance(candidate,new_main_rod) - 2*R0) > 1e-6 % adjust tolerence?
                        %                             candidate = CircularRod(contact_foots_j(i,:),triad,rc,R0,L0,phi_min_p,phi_max_p);
                        %                             deficit = 2*R0 - get_distance(candidate,new_main_rod);
                        %                             contact_foots_j(i,:) = contact_foots_j(i,:) + 1;
                        %                             candidate = CircularRod(contact_foots_j(i,:),triad,rc,R0,L0,phi_min_p,phi_max_p);
                        %                         end
                        
                        if get_distance(candidate,new_main_rod) < 2*R0;
                            ; % go on
                        elseif is_overlap(candidate,contacting_rods)
                            ; % go on
                        else
                            contacting_rods(i) = candidate;
                            break;
                        end
                        inneriter = inneriter + 1;
                    end
                    
                    if inneriter == max_inneriter
                        break;
                    end
                end
                
                outeriter = outeriter + 1;
                ind_nonempty = ~isempty(contacting_rods);
                if all(ind_nonempty)
                    break;
                end
                disp(sprintf('Perturbing contacting rods...\t Iteration: %d',outeriter));
            end
            
            found1_or_not0 = all(~isempty(contacting_rods));
            
            new_FRN.contacting_rods = contacting_rods;
            new_FRN.all_rods = [new_main_rod,contacting_rods];
            
            i = ones(num_contacts,1);
            j = (2:num_contacts+1)';
            
            contact_table = table(i,j,contact_points,contact_foots,contact_foots_j,...
                phi_contacts,theta_contacts);
            contact_table = renamevars(contact_table,...
                contact_table.Properties.VariableNames, ...
                ["i","j","p_i","s_i","s_j","phi_c","theta_c"]);
            
            new_FRN.contact_table = contact_table;
        end
        
        function ctable = get_contact_table(obj)
            N = numel(obj.contacting_rods);
            i = zeros(N,1);
            j = zeros(N,1);
            p_i = zeros(N,3);
            s_i = zeros(N,3);
            s_j = zeros(N,3);
            Q_i = zeros(N,9);
            Q_j = zeros(N,9);
            for i_c = 1:N
                [d_min,x1,x2] = get_distance(obj.main_rod,obj.contacting_rods(i_c));
                foot_m = obj.main_rod.query_points(x1);
                foot_c = obj.contacting_rods(i_c).query_points(x2);
                
                %                 close all;
                %                 plot(obj.main_rod,'r-');hold on;
                %                 plot(obj.contacting_rods(i_c),'b-');
                %                 plot3v(foot_m,'ro');
                %                 plot3v(foot_c,'bo');
                
                assert( abs(d_min - norm(foot_m - foot_c)) < 1e-12 );
                
                i(i_c) = 1;
                j(i_c) = i_c+1;
                p_i(i_c,:) = (foot_m + foot_c)/2;
                s_i(i_c,:) = foot_m;
                s_j(i_c,:) = foot_c;
                
                dir_m = obj.main_rod.query_directors(x1);
                u1_m = (foot_m - foot_c)/norm(foot_m - foot_c);
                u2_m = cross(dir_m(7:9),u1_m);u2_m = u2_m/norm(u2_m);
                u3_m = cross(u1_m,u2_m);u3_m = u3_m/norm(u3_m);
                Q_i(i_c,:) = [u1_m,u2_m,u3_m];
                
                dir_c = obj.contacting_rods(i_c).query_directors(x2);
                u1_c = (foot_c - foot_m)/norm(foot_m - foot_c);
                u2_c = cross(dir_c(7:9),u1_c);u2_c = u2_c/norm(u2_c);
                u3_c = cross(u1_c,u2_c);u3_c = u3_c/norm(u3_c);
                Q_j(i_c,:) = [u1_c,u2_c,u3_c];
            end
            ctable = table(i,j,p_i,s_i,s_j,Q_i,Q_j);
        end
        
        function [image_stacks,left_corner,right_corner] = export_image_stacks(obj_array)
            N = numel(obj_array);
            left_corners = zeros(N,3);
            right_corners = zeros(N,3);
            
            voxellist_to_export = cell(N,1);
            for i = 1:N
                rod_radius = obj_array(i).all_rods(1).rod_radius;
                voxel_list = obj_array(i).all_rods.get_voxel_list(1000);
                fat_rods = cellfun(@(x) sphere_sweep(x,rod_radius*1.2),voxel_list,'uniformoutput',false);
                voxellist_to_export{i} = fat_rods;
                voxels = round(vertcat(fat_rods{:}));
                left_corners(i,:) = min(voxels);
                right_corners(i,:) = max(voxels);
            end
            
            left_corner = min(left_corners);
            right_corner = max(right_corners);
            
            left_corner = left_corner - [10,10,10];
            right_corner = right_corner + [20,20,20];
            image_size = right_corner - left_corner;
            
            image_stacks(N).stack = zeros(image_size,'logical');
            
            for i = 1:N
                voxel_cell = cellfun(@(x) unique(round(x)-left_corner,'rows','stable'),voxellist_to_export{i},'uniformoutput',false);
                voxel_array = vertcat(voxel_cell{:});
                stack = zeros(image_size,'logical');
                ind = sub2ind2(image_size,voxel_array);
                stack(ind) = 1;
                image_stacks(i).stack = stack;
            end
            
        end
        
        function tf = isempty(obj_array)
            N = numel(obj_array);
            tf = zeros(N,1,'logical');
            for i = 1:N
                tf(i) = isempty(obj_array(i).main_rod);
            end
        end
        
        function h = plot(obj_array,varargin)
            N = numel(obj_array);
            for i = 1:N
                h = plot(obj_array(i).all_rods,varargin{:});hold on;
            end
        end
        
        function [Q_i,Q_j] = foo(obj)
            num_contacts = size(obj.contact_table,1);
            R0 = obj.main_rod.rod_radius;
            phi_contacts = obj.contact_table.phi_c;
            theta_contacts = obj.contact_table.theta_c;
            
            contact_foots = obj.main_rod.query_points(phi_contacts);
            directors = obj.main_rod.query_directors(phi_contacts);
            contact_points = contact_foots ...
                + R0*cos(theta_contacts).*directors(:,1:3) ...
                + R0*sin(theta_contacts).*directors(:,4:6);
            contact_foots_j = zeros(size(contact_foots));
            Q_i = zeros(num_contacts,9);
            Q_j = zeros(num_contacts,9);
            for i = 1:num_contacts
                v_ij = contact_points(i,:) - contact_foots(i,:);
                contact_foots_j(i,:) = contact_foots(i,:) + 2*v_ij; % contact foot at rod j
                
                n_ij = -v_ij/norm(v_ij);
                b_ij = cross(directors(i,7:9),n_ij);
                t_ij = cross(n_ij,b_ij);
                
                n_ji = -n_ij;
                b_ji = cross(obj.contacting_rods(i).plane_triad(7:9),n_ji);
                t_ji = cross(n_ji,b_ji);
                
                Q_i(i,:) = [n_ij,b_ij,t_ij];
                Q_j(i,:) = [n_ji,b_ji,t_ji];
            end
            
        end
    end
end

function out = uniform_log_random(n,a,b)
out = 10^( a*rand(n,1) + b);
end

function tf = is_overlap(rod_i,rod_array)
tf = 0;
N = numel(rod_array);
Ri = rod_i.rod_radius;
for j = 1:N
    rod_j = rod_array(j);
    Rj = rod_j.rod_radius;
    d_ij = get_distance(rod_i,rod_j);
    if d_ij < (Ri + Rj)*1.2
        tf = 1;
        return;
    end
end
end


function sweep = sphere_sweep(centerline,r)
    % get voxel indices from inequalities
%     assert(mod(r,2) == 1)
    centerline = unique(round(centerline),'rows','stable');    

    ind = (1: (2*r+1) ) - (r+1);
    [X,Y,Z] = meshgrid(ind,ind,ind);
    I = X.^2 + Y.^2 + Z.^2 <= r^2;
    sphere_points = ind2sub2(size(I),find(I)) - (r+1);
    
    num_centerline_points = size(centerline,1);
    num_sphere_points = size(sphere_points,1);
    
    all_points = zeros(num_centerline_points*num_sphere_points,3);
%     num_sphere_points*(1-1)+num_sphere_points+1;
    for i = 1:num_centerline_points        
%         size( all_points(num_sphere_points*(i-1)+i:num_sphere_points*(i-1)+num_sphere_points,:) )
%         size( sphere_points)
        
        all_points(num_sphere_points*(i-1)+1:num_sphere_points*(i-1)+num_sphere_points,:)...
            = centerline(i,:) + sphere_points;
    end
    sweep = unique(all_points,'rows','stable');    
end