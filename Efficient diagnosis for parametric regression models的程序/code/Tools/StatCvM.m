function [Tn, A] = StatCvM(X, e)
%%%------Compute Cvm statistic------%%% 
% input:  X      n*p  sample of covariates
%         e      n*1  residuals
% output: Tn     1*1  test statistic
%         A      n*n  A(i,j) = \sum_{k=1}^n A_{ijk}
[n, p] = size(X);
C_p = pi^(p/2-1)/gamma(p/2+1); 
A = zeros(n,n); 
temp = cell(n,1);
for i = 1 : n
    temp{i,1} = repmat(X(i,:),n,1) - X;
end
if p > 1
    for i = 1 : n
        A(i,i) = C_p * (n+1) * pi;
        for j = i+1:n        
            nume = sum(temp{i,1} .* temp{j,1},2);
            deno = sqrt(sum(temp{i,1}.^2,2)) .* sqrt(sum(temp{j,1}.^2,2));
            nume(deno==0) = 1; deno(deno==0) = 1; % 0/0=1 arcos(1)=0
            A(i,j) = C_p * sum( abs( pi-acos(nume./deno) ) );    
            A(j,i) = A(i,j);
        end
    end
elseif p==1
    for i = 1 : n
        for j = i:n        
            A(i,j) = sum((temp{i,1} <= 0) .* (temp{j,1} <= 0));
            A(j,i) = A(i,j);
        end
    end
end
Tn = e' * A * e / n^2; 
end