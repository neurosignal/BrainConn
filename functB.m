%% Assignment3: NBEHBC course (Pre-assignment:4)
%  Student: Amit K Jaiswal 
%  Date: 22-05-2017
% The function receives as input two vectors and a value corresponding to a
% number of permutations. The function should then, at each permutation, 
% swap a subset of values between the two vectors and compute the mean difference.
% The function should return the list of surrogate values. It also plots the 
% distribution of the surrogates (histogram) as well as where the unpermuted 
% difference of the means is (as a red vertical line). 
function [sv_dist, p_val] = functB(vector1,vector2, nperm)

if ~isequal(length(vector1),length(vector2))
   display('Input vector dimension mismatch');
   return
end

for ii=1:nperm;
    rand_series1=randi([0,1],1,length(vector1));
    rand_series1=rand_series1*2-1;
    rand_series2=-rand_series1;    
    sv_dist(ii)= mean(rand_series1.*vector1)+ mean(rand_series2.*vector2);
end

p_val = 1-sum(sv_dist<mean(vector2-vector1))/nperm;
display(['p_value=' num2str(p_val)]);
    
figure('color', [1 1 1])
histogram(sv_dist)
hold on
YLim = get(gca, 'ylim');
plot([mean(vector2-vector1)  mean(vector2-vector1)], YLim,'r','Linewidth', 2)
title('Mean difference distribution');
    
end

