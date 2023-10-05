function bw = binarymatrix(CC)

bw = zeros(CC.ImageSize,'logical');

for k = 1 : CC.NumObjects
    bw(CC.PixelIdxList{k}) = 1;
end

end