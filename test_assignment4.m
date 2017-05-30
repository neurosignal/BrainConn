%% Assignment3: NBEHBC course (Pre-assignment:4)
%  Student: Amit K Jaiswal 
%  Date: 22-05-2017
%% Define input arguments
vector=rand(1,20)*25;
vector1=rand(1,20)*35;
vector2=rand(1,20)*45;
nperm = 100;

%% Test function A
[sv_dist1, p_val1] = functA(vector,nperm);

%% Test function B
[sv_dist2, p_val2] = functB(vector1,vector2,nperm);
