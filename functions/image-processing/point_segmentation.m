function bins = point_segmentation(rr,varargin)

if nargin == 2
    threshold = varargin{1};
elseif nargin == 1
    threshold = 1.74;
end

dmat = pdist2(rr,rr);
amat = dmat < threshold;
G = graph(amat);
bins = conncomp(G);

end