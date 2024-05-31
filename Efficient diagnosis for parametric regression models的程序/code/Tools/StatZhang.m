function [Zn] = StatZhang(l_X, e)
%%%------Compute statistic: Zhang(2015)------%%% 
% input:  l_X    n*1  weight function  
%         e      n*1  residuals
% output: Zn     1*1  test statistic
n = size(e, 1);
II = (repmat(e,1,n) <= repmat(e',n,1));
Zn = sum((II' * l_X ).^2) / n^2; 
end