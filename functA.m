%% Assignment3: NBEHBC course (Pre-assignment:4)
%  Student: Amit K Jaiswal 
%  Date: 21-05-2017
function [sv_dist] = functA(vector,nperm)
for a=1:nperm;
    rand_series=randi([0,1],1,length(vector));
    rand_series=rand_series*2-1;
    surrog_vector=rand_series.*vector;
    sv_dist(1,a)=mean(surrog_vector);
end

    p_val = 1 - sum(mean(sv_dist) < mean(vector)) / nperm;
    display(['p_value=' num2str(p_val)]);
    
    figure('color', [1 1 1])
    histogram(sv_dist);
    hold on
    YLim = get(gca, 'ylim');
    plot([mean(sv_dist) mean(sv_dist)], YLim,'r','linewidth', 2);
    title('Surrogate vector mean distribution');
    
end
%%%%%%%%%%%% END