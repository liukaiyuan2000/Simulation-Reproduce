function [power] = RejectRate(Tn, Tn_star, alpha)
%%%------Compute rejection rates------%%%
% input:  Tn       1*N    test statistics 
%         Tn_star  B*N    bootstrap test statistics
%         alpha    1*int  significant level 
% output: power    int*1  empirical rejection rates
B = size(Tn_star,1);
STn = sort(Tn_star); % Sort each column from smallest to largest? B*N 
power = mean(repmat(Tn,length(alpha),1) >= STn(round((1-alpha)*B),:), 2); % line average: length(alpha) * N
end