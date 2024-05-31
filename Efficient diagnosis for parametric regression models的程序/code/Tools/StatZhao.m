function [Jn] = StatZhao(X, e)
%%%------Compute statistic: Zhao(2018)------%%% 
% input:  X      n*p  sample of covariates
%         e      n*1  residuals
% output: Jn     1*1  test statistic
[n, p] = size(X);
h = 0.5 * n^(-1/8); 
K = ones(n,n); % productive kernel function
for i = 1 : p
    XX =  repmat(X(:,i),1,n) - (repmat(X(:,i),1,n))';
    K = K .* (3/4) * (1-(XX/h).^2) .* (abs(XX)<=h);
end      
Vn = e' * K * e - (e.^2)' * diag(K); % leave-one-out estimator
Qn = 2 * ((e.^2)' * (K.^2) * (e.^2) - (e.^4)' * diag(K.^2)); % leave-one-out estimator
Jn = sqrt(n/(n-1)) * Vn / max(sqrt(Qn),1e-6);
end