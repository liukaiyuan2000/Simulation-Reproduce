function [Tn] = StatCvMRand(X, theta, e)
%%%------Compute Cvm statistic: random approximation------%%% 
% input:  X      n*p  sample of covariates
%         theta  p*M  random project direction
%         e      n*B  residuals
% output: Tn     1*B  test statistic: each element is based on the corresponding column of e
[n, B] = size(e);
M = size(theta, 2);
Tn = zeros(M,B);
Xt = X*theta;
for m = 1:M
    II = (repmat(Xt(:,m),1,n) <= repmat(Xt(:,m)',n,1));
    Tn(m,:) = sum((II' * e ).^2) / n^2; 
end
Tn = mean(Tn);
end