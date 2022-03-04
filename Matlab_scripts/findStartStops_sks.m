function [ starts, stops ] = findStartStops(t, gap_thr)


% takes t (a list of samples above the noise threshold) and finds the
% differences between adjacent numbers. If the difference between two
% adjacent numbers is greater than the gap threshold, record this as the
% end of a song. 

% resort the t list and do it again finding the spots where songs begin


    stops = [t(diff(t)>gap_thr)' t(end)];
    t_r = sort(t,'descend');
    starts = [t_r(diff(t_r)<-gap_thr)' t(1)];
    starts = sort(starts);
    
end

