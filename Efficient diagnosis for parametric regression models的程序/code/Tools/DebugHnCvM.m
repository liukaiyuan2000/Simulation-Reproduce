function [power, EstiPara] = DebugHnCvM(n, C, setting, alpha, N, B, hn)
%%% Test procedure of this paper
%%% debug bandwidth h =  hn * n^(-1/3) for CvM 
%% ------------- Initialization ------------- %%
switch setting
    case 1   % Zhang's example 1-model 1: LM p=2 q=0 
        beta = [2;3];
    case 2   % Zhang's example 1-model 2: NLM p=2 q=0 
        beta = [1;2]; 
    case 3   % LM p=5 q=0
        beta = [1;1;1;1;1];
    case 4   % NLM p=3 q=2  
        beta = [1;1;-1;-1;-1];
    case 5   % LM p=6 q=4
        beta = [1;1;1;1;-1;-1;-1;-1;-1;-1];
end
p = length(beta); 
%---statistics----bootstrap statistics------estimator of beta-----adjusted variables-----%
Tn = zeros(1,N); Tn_star = zeros(B,N); Parabeta  = zeros(N,p); hX = zeros(n,p);   % CvM
%%
for tt = 1 : N
    %% ---------- Generate data ---------- %%
    [tildeX, tildeY, Z, U] = GeneData(n, C, beta, setting);
    %% ------------ Error correction ------------ %%
    UU = repmat(U,1,n) - (repmat(U,1,n))';
    %------Local linear------%
    h =  hn * std(U) * n^(-1/3); 
    Kh = 3/4 * (1-(UU/h).^2) .* (abs(UU)<=h) / h; % standardized Ep kernel function
    hY = CaliLocaLine(tildeY, UU, Kh); 
    for i = 1 : size(tildeX,2)
        hX(:,i) = CaliLocaLine(tildeX(:,i), UU, Kh);
    end
    if size(Z, 2) > 0
        hX(:,p-size(Z,2)+1:p) = Z;
    end
    %% ---- Estimator of beta: (N)LSE ---- %%       
    [Parabeta(tt,:), e] = EstiBeta(hX, hY, setting); % this paper
    %% ---------- Compute test statistics ---------- %% 
    [Tn(tt), A] = StatCvM(hX, e);                       % CvM 
    %% --------- Bootstrap resampling --------- %%
    r = unifrnd(0,1,n,B);
    V = (1+sqrt(5))/2 * (r<(5-sqrt(5))/10) + (1-sqrt(5))/2 * (r>=(5-sqrt(5))/10);
    e_star = BootResi(hX, hY, e, V, setting);          % bootstrap residuals   CvM 
    Tn_star(:,tt) = diag(e_star' * A * e_star / n^2);  % bootstrap statistics  CvM                     
end % tt in 1:N
%% ---------- Compute rejection rates ---------- %%
power(1,:) = RejectRate(Tn, Tn_star, alpha);                       % CvM 
EstiPara = [mean(Parabeta)-beta'; sum(bsxfun(@minus,Parabeta,beta').^2)/N; std(Parabeta)];  % [bias; mse; std]
end