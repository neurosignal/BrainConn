%% Assignment3: NBEHBC course (Pre-assignment:4)
%  Student: Amit K Jaiswal 
%  Date: 21-05-2017
% The function receives a vector and a value corresponding to the number of
% permutations. It returns a list of mean surrogate values and p-value. It 
% also plots surrogate vector distribution histogram as well as where the actual 
% mean is (as a red vertical line).
function [sv_dist, p_val] = functA(vector,nperm)
for ii=1:nperm;
    rand_series=randi([0,1],1,length(vector));
    rand_series=rand_series*2-1;
    surrog_vector=rand_series.*vector;
    sv_dist(ii)= mean(surrog_vector);
end

    p_val = 1 - sum(sv_dist < mean(vector)) / nperm;
    display(['p_value=' num2str(p_val)]);
    
    figure('color', [1 1 1])
    histogram(sv_dist);
    hold on
    YLim = get(gca, 'ylim');
    plot([mean(sv_dist) mean(sv_dist)], YLim,'r','linewidth', 2);
    title('Surrogate vector mean distribution');
    
end
%%%%%%%%%%%% END