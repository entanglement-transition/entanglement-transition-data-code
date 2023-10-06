function new_edge = clip_edge_y(edge,y1,y2)
    if size(edge,1) < 1
        new_edge = edge;
        return;
    end
    if edge(3) == edge(6)
        new_edge = edge;
        return;
    end

    % single edge
    if edge(3) < edge(6)
        u = edge(1:3);
        v = edge(4:6);
    else
        u = edge(4:6);
        v = edge(1:3);
    end
    z = (y1+y2)/2;
    w = (y2-y1)/2;
    
    interpolation_factor1 = (y1 - u(2))/(v(2) - u(2));
    interpolation_factor2 = (y2 - v(2))/(u(2) - v(2));

    interpolation_factor1 = clip(interpolation_factor1,0,1);
    interpolation_factor2 = clip(interpolation_factor2,0,1);

    u_new = u + interpolation_factor1 * (v - u);
    v_new = v + interpolation_factor2 * (u - v);
    
    new_edge = [u_new,v_new];
end

% function new_edge = clip_edge(edge,y1,z2)
%     if size(edge,1) < 1
%         new_edge = edge;
%         return;
%     end
%     % single edge
%     if edge(3) < edge(6)
%         u = edge(1:3);
%         v = edge(4:6);
%     else
%         u = edge(4:6);
%         v = edge(1:3);
%     end
%     z = (y1+z2)/2;
%     w = (z2-y1)/2;
%     height_difference = abs(u(3)-v(3));
%     
%     if u(3) == v(3)
%         new_edge = [u,v];
%         return
%     end
%     
%     if (abs(u(3) - z) < w) && (abs(v(3) - z) < w) % both are in the layer
%         u_new = u;
%         v_new = v;
%     elseif (abs(u(3) - z) < w) && (abs(v(3) - z) > w) % u in the layer
%         u_new = u;
%         x2 = (v(1) - u(1))*(z2-u(3))/(v(3)-u(3)) + u(1);
%         y2 = (v(2) - u(2))*(z2-u(3))/(v(3)-u(3)) + u(2);
%         v_new = [x2,y2 z2];
%     elseif (abs(u(3) - z) > w) && (abs(v(3) - z) < w) % v in the layer
%         x1 = -(u(1) - v(1))*(y1-v(3))/(u(3)-v(3)) + v(1);
%         y1 = -(u(2) - v(2))*(y1-v(3))/(u(3)-v(3)) + v(2);
%         u_new = [x1,y1 y1];
%         v_new = v;
%     else % both are out of the layer
%         x1 = -(u(1) - v(1))*(y1-v(3))/(u(3)-v(3)) + v(1);
%         y1 = -(u(2) - v(2))*(y1-v(3))/(u(3)-v(3)) + v(2);
%         u_new = [x1,y1 y1];
%         
%         x2 = (v(1) - u(1))*(z2-u(3))/(v(3)-u(3)) + u(1);
%         y2 = (v(2) - u(2))*(z2-u(3))/(v(3)-u(3)) + u(2);
%         v_new = [x2,y2 z2];
%     end
%     new_edge = [u_new,v_new];
% end