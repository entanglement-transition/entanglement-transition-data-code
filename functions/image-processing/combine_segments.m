function [rr_new,end_to_end] = combine_segments(rr1,rr2)

ep1_1 = rr1(1,:);
ep1_2 = rr1(end,:);

ep2_1 = rr2(1,:);
ep2_2 = rr2(end,:);

distances = [norm(ep1_2 - ep2_1),
    norm(ep1_2 - ep2_2),
    norm(ep1_1 - ep2_1),
    norm(ep1_1 - ep2_2)];

[~,I] = min(distances);

switch I
    case 1
        rr_new = [rr1; rr2];
        end_to_end = [rr1(end,:);rr2(1,:)];
    case 2
        rr_new = [rr1; flip(rr2,1)];
        end_to_end = [rr1(end,:);rr2(end,:)];
    case 3
        rr_new = [flip(rr1,1); rr2];
        end_to_end = [rr1(1,:);rr2(1,:)];
        
    case 4
        rr_new = [flip(rr1,1); flip(rr2,1)];
        end_to_end = [rr1(1,:);rr2(end,:)];
    otherwise
        error('cannot find proper min');
end

end