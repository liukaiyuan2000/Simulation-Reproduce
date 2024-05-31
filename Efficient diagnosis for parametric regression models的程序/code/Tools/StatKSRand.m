function [Kn] = StatKSRand(X, theta, e)
%%%------Compute KS statistic: random approximation, absolute sign is inside the integral------%%% 
% input:  X      n*p  sample of covariates
%         theta  p*M  random project direction
%         e      n*B  residuals
% output: Kn     1*B  test statistic: each element is based on the corresponding column of e
[n, B] = size(e);
M = size(theta, 2);
Kn = zeros(1,B);
Xtheta = X * theta; % n*M
Bv = zeros(n,M,B);
for j = 1 : n
    for k = 1 : M
        Bv(j,k,:) = sum(abs(e' * double(Xtheta <= Xtheta(j,k))), 2);
    end
end
Kn(1,:) = max(max(Bv)) / (M * sqrt(n));
end



