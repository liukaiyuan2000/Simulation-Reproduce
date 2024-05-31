function [hatbeta, e] = EstiBeta(X, Y, setting)
%%%------Estimate beta: (N)LSE------%%%  
% input:  X        n*p  sample of covariates
%         Y        n*1  sample of response
%         setting  1*1  simulation settings
% output: hatbeta  p*1  estimated parameters
%         e        n*1  residuals
p = size(X,2);
if setting == 2 % NLSE
    options = optimoptions('lsqcurvefit','Display','off'); 
    hatbeta = lsqcurvefit(@(x,xdata)(x(1) + xdata(:,1).*((1+xdata(:,2)).^(x(2)))), zeros(p,1), X, Y, [], [], options);
    e = Y - hatbeta(1) - X(:,1).*((1+X(:,2)).^(hatbeta(2))); % residual n*1
elseif setting == 4 % NLSE
    options = optimoptions('lsqcurvefit','Display','off'); 
    hatbeta = lsqcurvefit(@(x,xdata)(exp(xdata*x)), zeros(p,1), X, Y, [], [], options);
    e = Y - exp(X*hatbeta); 
else % LSE
    hatbeta = (X'*X+1e-6)\X'*Y; 
    e = Y - X*hatbeta; 
end
end