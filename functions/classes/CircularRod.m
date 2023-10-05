classdef CircularRod
    properties
        rod_radius
        rod_length
        plane_triad
        plane_origin
        
        radius_curvature % radius of curvature
        phi_min
        phi_max
    end
    
    methods
        function obj = CircularRod(plane_origin,...
                plane_triad,... % centerline pass through this point
                radius_curvature,...
                rod_radius,...
                rod_length,phi_min,phi_max)
            
            if nargin == 0;
                ;
            else nargin > 0;
                obj.plane_triad = plane_triad;
                obj.plane_origin = plane_origin;
                obj.rod_length = rod_length;
                obj.rod_radius = rod_radius;
                obj.radius_curvature = radius_curvature;
                obj.phi_min = phi_min;
                obj.phi_max = phi_max;
                
                assert( abs( (obj.phi_max - obj.phi_min)*obj.radius_curvature - obj.rod_length )  < 1e-12 );
                assert( obj.phi_max >= 0 );
                assert( obj.phi_min <= 0 );
            end
        end
        
        
        function random_perturb(obj,v,Q,fixed_point)
            % v: translation vector
            % Q: rotation
            
        end
        
        function [spatial_points,philist] = get_discrete_points(obj,N)
            philist = linspace(obj.phi_min,obj.phi_max,N)';
            xpt = +obj.radius_curvature*sin(philist);
            ypt = -obj.radius_curvature*cos(philist) + obj.radius_curvature;
            
            u1 = obj.plane_triad(1:3);
            u3 = obj.plane_triad(7:9);
            spatial_points = xpt.*u3 + ypt.*u1;
            spatial_points = spatial_points + obj.plane_origin;
        end
        
        function directors = query_directors(obj,philist)
            assert( min(philist) >= obj.phi_min );
            assert( max(philist) <= obj.phi_max );
            
            xv = obj.radius_curvature*cos(philist);
            yv = obj.radius_curvature*sin(philist);
            
            u1 = obj.plane_triad(1:3); % tangent vector at the origin
            u2 = obj.plane_triad(4:6); % plane normal
            u3 = obj.plane_triad(7:9); % normal to curve in the osculating plane
            
            d3 = xv.*u3 + yv.*u1; % tangent at query points
            d3 = d3./rwnorm(d3);
            
            d2 = repmat(u2, [numel(philist),1]);
            
            d1 = cross(d2,d3,2);
            d1 = d1./rwnorm(d1);
            
            directors = [d1,d2,d3]; % 1,9 row vector
            
        end
        
        function [tangents,philist] = get_discrete_tangents(obj,N)
            philist = linspace(obj.phi_min,obj.phi_max,N)';
            
            xv = obj.radius_curvature*cos(philist);
            yv = obj.radius_curvature*sin(philist);
            
            u1 = obj.plane_triad(1:3);
            u3 = obj.plane_triad(7:9);
            tangents = xv.*u3 + yv.*u1;
        end
        
        function spatial_points = query_points(obj,philist)
            assert( min(philist) >= obj.phi_min );
            assert( max(philist) <= obj.phi_max );
            
            xpt = obj.radius_curvature*sin(philist);
            ypt = -obj.radius_curvature*cos(philist) + obj.radius_curvature;
            
            u1 = obj.plane_triad(1:3);
            u3 = obj.plane_triad(7:9);
            spatial_points = xpt.*u3 + ypt.*u1;
            spatial_points = spatial_points + obj.plane_origin;
            
        end
        
        function [d_min,phi_opt1,phi_opt2] = get_distance(obj1,obj2)
            if isempty(obj1) | isempty(obj2)
                d_min = NaN;
                phi_opt1 = NaN;
                phi_opt2 = NaN;
                return;
            end
            rc1 = obj1.radius_curvature;
            Q1 = obj1.plane_triad;
            p1 = obj1.plane_origin;
            phi1_l = obj1.phi_min;
            phi1_u = obj1.phi_max;
            
            rc2 = obj2.radius_curvature;
            Q2 = obj2.plane_triad;
            p2 = obj2.plane_origin;
            phi2_l = obj2.phi_min;
            phi2_u = obj2.phi_max;
            
            dist_fun = @(th) norm(p1+(rc1*sin(th(1)).*Q1(7:9) - (rc1*cos(th(1))-rc1).*Q1(1:3))...
                -( p2+(rc2*sin(th(2)).*Q2(7:9)-(rc2*cos(th(2))-rc2).*Q2(1:3)) ) );
            
            %             A = [[1,0];[-1,0];[0,1];[0,-1]];
            %             b = [phi1_u;-phi1_l;phi2_u;-phi2_l];
            
            lb = [phi1_l,phi2_l]; ub = [phi1_u,phi2_u];
            [pts1,phi1] = obj1.get_discrete_points(60);
            [pts2,phi2] = obj2.get_discrete_points(60);
            
            dmat = pdist2(pts1,pts2);
            [quick_min,i_min] = min(dmat,[],'all','linear');
            [i1,i2] = ind2sub(size(dmat),i_min);
            x0 = [phi1(i1),phi2(i2)];
            options = optimoptions('fmincon','Display','off','Algorithm','sqp','ConstraintTolerance',1e-12,'StepTolerance',1e-12);
            xopt = fmincon(dist_fun,x0,[],[],[],[],lb,ub,[],options);
            
            d_min = dist_fun(xopt);
            phi_opt1 = xopt(1);
            phi_opt2 = xopt(2);
        end
        
        function [d_min,phi_opt] = get_distance_from_point(obj,pt)
            rc1 = obj.radius_curvature;
            Q1 = obj.plane_triad;
            p1 = obj.plane_origin;
            phi1_l = obj.phi_min;
            phi1_u = obj.phi_max;
            lb = [phi1_l]; ub = [phi1_u];

            dist_fun = @(th) norm(p1+(rc1*sin(th(1)).*Q1(7:9) - (rc1*cos(th(1))-rc1).*Q1(1:3))...
                -( pt ) );

            [pts1,phi1] = obj.get_discrete_points(60);
            [quick_min,i_min] = min(norm(pts1 - pt),[],'all','linear');
            
            x0 = phi1(i_min);
            options = optimoptions('fmincon','Display','off','Algorithm','sqp','ConstraintTolerance',1e-12,'StepTolerance',1e-12);
            xopt = fmincon(dist_fun,x0,[],[],[],[],lb,ub,[],options);
            
            d_min = dist_fun(xopt);
            phi_opt = x0;
        end
        
        function voxel_list = get_voxel_list(obj_array,num_points)
            N = numel(obj_array);
            voxel_list = cell(N,1);
            for i = 1:N
                voxel_list{i} = obj_array(i).get_discrete_points(num_points);
            end
        end
        
        function h = plot(obj_array,varargin)
            N = numel(obj_array);
            for i = 1:N
                pts = obj_array(i).get_discrete_points(100);
                h = plot3v(pts,varargin{:});hold on;
            end
        end
        
        function tf = isempty(obj_array)
            N = numel(obj_array);
            tf = zeros(N,1,'logical');
            for i = 1:N
                tf(i) = isempty(obj_array(i).rod_radius);
            end
        end
    end
    
    methods (Static)
        
    end
end