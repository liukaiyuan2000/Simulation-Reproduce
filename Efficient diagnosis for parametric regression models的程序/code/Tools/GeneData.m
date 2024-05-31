function [tildeX, tildeY, Z, U] = GeneData(n, C, beta, setting)
%%%------Generate data------%%%
% input:  n        1*1  sample size 
%         C        1*1  distance between H0 and H1
%         beta     p*1  parameters
%         setting  1*1  simulation settings
% output: tildeX   n*   observed/distorted covariates
%         tildeY   n*1  observed/distorted response
%         Z        n*   observed covariates
%         U        n*1  confunding variable
epsilon = normrnd(0,1,n,1); % random error
U = unifrnd(0,1,n,1); % confounding variable
switch setting
    case 1  % Zhang's example 1-model 1: LM p=2 q=0
        X = unifrnd(1,2,n,length(beta)); % real covariates: n*p  p+q=2 LM   
        tildeX = [(1 + 0.3*cos(2*pi*U)),(1 + 0.2*(U.^2-1/3))] .* X; % distorted covariates: n*p
        Z = [];
        Y = X*beta + 2*C*exp(0.5*X(:,2)) + 0.15 * epsilon; % real response: n*1
        tildeY = (1 + 0.2*cos(2*pi*U)) .* Y; % distorted response: n*1        
    case 2  % Zhang's example 1-model 2: NLM p=2 q=0 
        X = unifrnd(1,2,n,length(beta)); % p+q=2 NLM     
        tildeX = [(1 + 0.3*cos(2*pi*U)),(1 + 0.2*(U.^2-1/3))] .* X; 
        Z = [];
        Y = beta(1) + X(:,1).*((1+X(:,2)).^(beta(2))) + C*exp(0.5*X(:,2)) + 0.15 * epsilon; 
        tildeY = (1 + 0.2*cos(2*pi*U)) .* Y;        
    case 3  % LM p=5 q=0
        X = unifrnd(1,2,n,length(beta)); % p+q=5 LM
        tildeX = [(1 + 0.3*cos(2*pi*U)),(1 + 0.2*(U.^2-1/3)),(U + 1/2),(1 + 0.2*(U.^2-1/3)),(U.^2 + 2/3)] .* X;  
        Z = [];
        Y = X*beta + 2*C*exp(0.5*X(:,2)) + 0.15 * epsilon;
        tildeY = (1 + 0.2*cos(2*pi*U)) .* Y;        
    case 4 % NLM p=3 q=2        
        X = unifrnd(1,2,n,length(beta)); % p+q=5 NLM 
        tildeX = [(1 + 0.3*cos(2*pi*U)),(1 + 0.2*(U.^2-1/3)),(U.^2 + 2/3)] .* X(:,1:3);
        Z = X(:,size(tildeX,2)+1:end);
        Y = exp(X*beta) + C*(X(:,1:3)*beta(1:3)) + 0.15 * epsilon; 
        tildeY = (1 + 0.2*cos(2*pi*U)) .* Y;
    case 5 % LM p=6 q=4
        X = unifrnd(1,2,n,length(beta)); % p+q=10 LM  
        tildeX = [(1 + 0.3*cos(2*pi*U)),(1 + 0.3*cos(2*pi*U)),(1 + 0.3*cos(2*pi*U))...
            (1 + 0.2*(U.^2-1/3)),(1 + 0.2*(U.^2-1/3)),(1 + 0.2*(U.^2-1/3))] .* X(:,1:6);
        Z = X(:,size(tildeX,2)+1:end);
        Y = X*beta + C*0.1*exp(X(:,1)+X(:,2)) + 0.15 * epsilon;  
        tildeY = (1 + 0.2*cos(2*pi*U)) .* Y;        
end
end

% zhao2018 model 1
% X = normrnd(0,1,n,length(beta));     
% tildeX = [(U + 1/2),(1 + 0.2*(U.^2-1/3))] .* X; 
% Y = X*beta + C*exp(0.6*X(:,2))/2 + 0.5 * epsilon;
% tildeY = (U.^2 + 2/3) .* Y; 
% zhao 2018 model 2
% X = normrnd(0,1,n,length(beta));     
% tildeX = [(U + 1/2),(1 + 0.2*(U.^2-1/3))] .* X; 
% Y = 1.5 * sin(X*beta) + C*exp(0.6*X(:,2))/2 + 0.5 * epsilon;
% tildeY = (U.^2 + 2/3) .* Y; 