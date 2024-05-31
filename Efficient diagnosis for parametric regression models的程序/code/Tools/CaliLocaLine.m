function [hatY] = CaliLocaLine(Y, UU, Kh)
%%%------Error correction: local linear------%%% 
% input:  Y        n*1  sample of distorted variable
%         UU       n*n  based on the confunding variable U
%         Kh       n*n  based on the scaled kernel function of UU
% output: hatY     n*1  adjusted variable
n = size(UU,1);
s0 = sum(Kh);        % 1*n: n[s0(U1,h),s0(U2,h),...,s0(Un,h)]
s1 = sum(UU.*Kh);    % 1*n: n[s1(U1,h),s1(U2,h),...,s1(Un,h)] 
s2 = sum(UU.^2.*Kh); % 1*n: n[s2(U1,h),s2(U2,h),...,s2(Un,h)] 
nume = sum((repmat(s2,n,1)-repmat(s1,n,1) .* UU) .* Kh .* repmat(Y,1,n));
deno = s0.*s2-s1.^2;
nume(deno==0) = 1; deno(deno==0) = 1; % 0/0=1
EY  = (nume ./ deno)'; % local linear kernal estimator: E(Y|U)
hatY = mean(Y) * Y ./ EY;  % adjusted variables 
end
% t0 = Kh * Y;
% t1 = UU' .* Kh * Y;
