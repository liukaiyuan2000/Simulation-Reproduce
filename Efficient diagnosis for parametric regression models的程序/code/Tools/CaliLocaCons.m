function [hatY] = CaliLocaCons(Y, Kh, mode)
%%%------Error correction: local constant------%%% 
% input:  Y        n*1  sample of distorted variable
%         Kh       n*n  based on the scaled kernel function
%         mode     1*1  0: Y; 1: abs(Y)
% output: hatY     n*1  adjusted variable
if mode == 0
    EY = Kh * Y ./ sum(Kh,2);  % local constant kernal estimator: E[ Y | U ]
    Y(EY==0) = 1; EY(EY==0) = 1; % 0/0=1
    hatY = mean(Y) * Y ./ EY;  % adjusted variables 
else 
    EY = Kh * abs(Y) ./ sum(Kh)';   % local constant kernal estimator: E[ |Y| | U ]
    Y(EY==0) = 1; EY(EY==0) = 1; 
    hatY = mean(abs(Y)) * Y ./ EY;  
end