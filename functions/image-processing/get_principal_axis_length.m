function stats = get_principal_axis_length(rr)
    if isempty(rr)
        stats = [];
        return;
    end

    if size(rr,1) == 1
        stats.Centroid = rr;
        stats.PrincipalAxisLength = [0 0 0];
        stats.Orientation = [0 0 0];
        stats.EigenValues = [0 0 0];
        stats.EigenVectors = [0 0 0];
        return;
    end
    
%     centroid = stats(k).Centroid-stats(k).BoundingBox(1:3)+0.5;
    centroid = mean(rr,1);
    centered = (rr - centroid);
    mu000 = sum(centered(:).^0);
    mu200 = sum( centered(:,1).^2 ) / mu000 + 1/12;
    mu020 = sum( centered(:,2).^2 ) / mu000 + 1/12;
    mu002 = sum( centered(:,3).^2 ) / mu000 + 1/12;
    mu110 = sum( centered(:,1).*centered(:,2) ) / mu000;
    mu011 = sum( centered(:,2).*centered(:,3) ) / mu000;
    mu101 = sum( centered(:,3).*centered(:,1) ) / mu000;
    
    numPoints = size(rr,1);
    covMat = [mu200 mu110 mu101; ...
              mu110 mu020 mu011; ...
              mu101 mu011 mu002]./numPoints;
    
    [U,S] = svd(covMat);
    [S,ind] = sort(diag(S), 'descend');
    
    U = U(:,ind);
    % Update U so that the first axis points to positive x
    % direction and make sure that the rotation matrix determinant
    % is positive
    if U(1,1) < 0
        U = -U;
        U(:,3) = -U(:,3);
    end
    
    [V,D] = eig(covMat);
    [D,ind] = sort(diag(D), 'descend');
    
    V = V(:,ind);           
    
    stats.Centroid = centroid;
    stats.PrincipalAxisLength = [4*sqrt(S(1)*numPoints) 4*sqrt(S(2)*numPoints) 4*sqrt(S(3)*numPoints)];
    stats.Orientation = rotm2euler(U);
    stats.EigenValues = D*numPoints;
    stats.EigenVectors = V;
end

function centralMoments = calculateCentralMoments(im,centroid,i,j,k)

[r,c,p] = size(im);
centralMoments = ((1:r)-centroid(2))'.^i * ((1:c)-centroid(1)).^j;
z = reshape(((1:p)-centroid(3)).^k,[1 1 p]);
centralMoments = centralMoments.*z.*im;
centralMoments = sum(centralMoments(:));

end

function eulerAngles = rotm2euler(rotm)
%ROTM2EULER Convert rotation matrix to Euler angles
%
%   eulerAngles = rotm2euler(rotm) converts 3x3 3D rotation matrix to Euler
%   angles
%
%   Reference:
%   ---------
%   
%   Ken Shoemake, Graphics Gems IV, Edited by Paul S. Heckbert,
%   Morgan Kaufmann, 1994, Pg 222-229.

% Scale factor to convert radians to degrees
k = 180 / pi;

cy = hypot(rotm(1, 1), rotm(2, 1));

if cy > 16*eps(class(rotm))
    psi     = k * atan2( rotm(3, 2), rotm(3, 3));
    theta   = k * atan2(-rotm(3, 1), cy);
    phi     = k * atan2( rotm(2, 1), rotm(1, 1));
else
    psi     = k * atan2(-rotm(2, 3), rotm(2, 2));
    theta   = k * atan2(-rotm(3, 1), cy);
    phi     = 0;                    
end

eulerAngles = [phi, theta, psi];
end