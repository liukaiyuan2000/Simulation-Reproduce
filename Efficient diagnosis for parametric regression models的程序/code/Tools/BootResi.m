function [e_star] = BootResi(X, Y, e, V, setting)
%%%------Compute bootstrap residuals------%%%
% input:  X        n*p  sample of covariates
%         Y        n*1  sample of response
%         e        n*1  residuals
%         V        n*B  random number
%         setting  1*1  simulation settings
% output: e_star   n*B  bootstrap residuals         
[n, B] = size(V);
Y_star = repmat(Y-e,1,B) + repmat(e,1,B) .* V;  % n*B  
e_star = zeros(n,B);
for l = 1 : B            
    [~, e_star(:,l)] = EstiBeta(X, Y_star(:,l), setting);                 
end
end    