function plot_edges(edges,varargin)
N = size(edges,1);
for i = 1:N
    r1 = edges(i,1:3);
    r2 = edges(i,4:6);
    plot3v([r1;r2],varargin{:});hold on;
end
end