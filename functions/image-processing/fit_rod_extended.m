function out = fit_rod_extended(rr,varargin)
    % outputs: cen, pts, r, u, philist, philist2, xm, ym, len, ref, err
    % inputs: rr, DEBUG_FLAG
    % rr is a 3xN matrix of points
    % DEBUG_FLAG is a boolean flag to turn on debugging plots
    if nargin > 1
        DEBUG_FLAG = varargin{1};
    else
        DEBUG_FLAG = false;
    end
    
    n=size(rr,2); % number of data points
    if isempty(rr)
        out.cen= [];
        out.pts = [];
        out.r = NaN;
        out.u = NaN;
        out.philist = [];
        out.philist2 = [];
        out.xyr = [];
        out.err = [];
        out.len = 0;
        out.ref = rr;
        return
    end
    
    if n == 1
        out.cen= rr';
        out.pts = rr';
        out.r = NaN;
        temp = rr/norm(rr);
        out.u = [temp temp temp];
        out.philist = [];
        out.philist2 = [];
        out.xyr = [];
        out.err = [];
        out.len = 0;
        out.ref = rr;
        return
    end
    
    if n == 2
        out.cen= mean(rr',1);
        out.pts = linspacev(rr(:,1)',rr(:,end)',1000);
        out.r = Inf;
        ori = (rr(:,2) - rr(:,1))/norm((rr(:,2) - rr(:,1)));
        out.u = [ori ori ori];
        out.philist = [];
        out.philist2 = [];
        out.xyr = [];
        out.err = [];
        out.len = norm(rr(:,2) - rr(:,1));
        out.ref = rr;
        return
    end
    
    % utilizing 'frenet_robust' code
    % https://www.mathworks.com/matlabcentral/fileexchange/47885-frenet_robust-zip
    
    % fit coalescing plane ------------------------------------------------
    centroid=mean(rr,2);
    rrwin_centered=rr-repmat(centroid,1,size(rr,2));
    % least-square fit of coalescing plane using singular-value-decomposition
    [U,S,~]=svd(rrwin_centered);
    % sort singular values
    [~,ind]=sort(diag(S),'descend');
    S=S(ind,ind);
    U=U(:,ind);
    S=S/sum(diag(S)); % normalize singular values
    u1=U(:,1); % first unit vector spanning coalescing plane: agrees with tt0
    u2=U(:,2); % second unit vectors spanning coalescing plane: estimate for nn
    u3=U(:,3); % normal vector of coalescing plane
    orientation = U(:,1)*sign(sum(U(:,1).*(rrwin_centered(:,end) - rrwin_centered(:,1)),1));
    slist = sum((rr - centroid).*orientation ,1);
    dlist = cwnorm(rr - (centroid + slist.*orientation));
    [~,i_sort] = sort(slist);
    rrwin_centered = rrwin_centered(:,i_sort);
    
    if DEBUG_FLAG
        close all;
        plot(rrwin_centered'*u1, rrwin_centered'*u2,'o-');axis equal;
        figure;
        plot3v(rrwin_centered');axis equal;
    end
    
    N = 1000;
    if max(dlist) < 1e-10 % if all points are on a line (need to make tolerance passed in?)
        slist = sum((rr-centroid).*u1,1);
        s1 = min(slist);
        s2 = max(slist);
    
        r1 = centroid + s1*u1;
        r2 = centroid + s2*u1;
        out.cen= centroid;
        out.pts = linspacev(r1',r2',1000);
        out.r = Inf;
        out.u = [u1 u2 u3];
        out.philist = slist;
        out.philist2 = linspace(min(s1),min(s2),1000);
        out.xm = 0;
        out.ym = 0;
        out.len = s2 - s1;
    
        best_estimation = centroid + u1.*slist;
        out.ref = best_estimation;
    
        err = (rr' - best_estimation').^2;
        out.err = (rwnorm(err));
        return
    end
    
    %r0 < 1e7 %& n > 100 % optimal value to check line?
    % fit circle ----------------------------------------------------------
    param=CircleFitByTaubin([rrwin_centered'*u1 rrwin_centered'*u2]);
    xm=param(1); ym=param(2); % center of coalescing circle in (u1,u2)-plane
    RR=centroid+xm*u1+ym*u2; % center of coalescing circle in 3d-space
    r0=abs(param(3)); % radius of circle
    
    if isnan(r0) | isinf(r0)
        out.cen= centroid;
        out.r = Inf;
        out.u = [u1,u2,u3];
        out.pts = rr';
        out.philist = NaN;
        out.philist2 = NaN;
        out.xm = NaN;
        out.ym = NaN;
        out.len = sum(rwnorm(diff(rr')));
        out.ref = rr';
        out.err = Inf;
        return;
        
    end
    
    
    %         if dot( rr(:,end) - rr(:,1),u1 ) < 0
    %             u1 = -u1;
    %         end
    %         if dot(RR - centroid,u2) < 0
    %             u2 = -u2;
    %         end
    
    u3 = cross(u1,u2);
    %         dot(u1,u3)
    assert( abs( dot(u1,u3) ) < 1e-5 );
    
    %         philist=sort(unwrap(atan2(rrwin_centered'*u2-ym,rrwin_centered'*u1-xm)));
    philist = (unwrap(atan2(rrwin_centered'*u2-ym,rrwin_centered'*u1-xm)));
    
    delta_th = (max(philist) - min(philist))/(N/4);
    min_th = min(philist) - delta_th*(N/4);
    max_th = max(philist) + delta_th*(N/4);
    
%     min_th = min(philist);
%     max_th = max(philist);
    
    th = linspace(min_th,max_th,N);
    xpt = xm+r0*cos(th);
    ypt = ym+r0*sin(th);
    
    xpt2 = xm+r0*cos(philist)';
    ypt2 = ym+r0*sin(philist)';    
    best_estimation = (centroid + xpt2.*u1 + ypt2.*u2);
    err = (rr - best_estimation).^2;
    
    circle_points3D = centroid + xpt.*u1 + ypt.*u2;    
    
    out.cen= centroid;
    out.r = r0;
    out.u = [u1,u2,u3];
    out.pts = circle_points3D';
    out.philist = philist;
    out.philist2 = th;
    out.xm = xm;
    out.ym = ym;
    out.len = r0*(max_th - min_th);
    out.ref = best_estimation;
    out.err = (rwnorm(err'));
    
    if norm(rr(:,1)' - out.pts(1,:)) > norm(rr(:,1)' - out.pts(end,:))
        out.pts = flipud(out.pts);
    end
    
    % visual check --------------------------------------------------------
    if DEBUG_FLAG
        %         figure
        %         plot(rrwin_centered'*u1, rrwin_centered'*u2,'o','linewidth',1);
        %         hold on
        %         plot(xpt,ypt,'linewidth',2);
        %         axis equal
        %             my_plot3(circle_points3D','-','linewidth',2)
        my_plot3(best_estimation','o-','linewidth',2);
        hold on
        my_plot3(rr','o');
        axis equal
    end;
    
    end
    