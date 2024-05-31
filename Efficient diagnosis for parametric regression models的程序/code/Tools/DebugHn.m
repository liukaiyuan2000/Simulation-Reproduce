function [power, EstiPara] = DebugHn(n, C, setting, alpha, N, B, M, hn)
%%% Test procedure of this paper
%%% debug bandwidth h =  hn * n^(-1/3) for all our test statistics
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
Tn  = zeros(1,N); Tn_star  = zeros(B,N); Parabeta   = zeros(N,p); hX = zeros(n,p);   % CvM
TnU = zeros(1,N); TnU_star = zeros(B,N);                                             % CvM: uniform approximation
TnN = zeros(1,N); TnN_star = zeros(B,N);                                             % CvM: gaussian approximation
KnU = zeros(1,N); KnU_star = zeros(B,N);                                             % KS: uniform approximation
KnN = zeros(1,N); KnN_star = zeros(B,N);                                             % KS: gaussian approximation
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
    theta1 = unifrnd(-1,1,p,M); 
    theta2 = mvnrnd(zeros(M,1),eye(M),p); 
    [Tn(tt), A] = StatCvM(hX, e);                       % CvM 
    TnU(tt) = StatCvMRand(hX, theta1, e);               % CvM: uniform approximation
    TnN(tt) = StatCvMRand(hX, theta2, e);               % CvM: gaussian approximation
    KnU(tt) = StatKSRand(hX, theta1, e);                % KS: uniform approximation
    KnN(tt) = StatKSRand(hX, theta2, e);                % KS: gaussian approximation
    %% --------- Bootstrap resampling --------- %%
    r = unifrnd(0,1,n,B);
    V = (1+sqrt(5))/2 * (r<(5-sqrt(5))/10) + (1-sqrt(5))/2 * (r>=(5-sqrt(5))/10);
    e_star = BootResi(hX, hY, e, V, setting);           % bootstrap residuals   CvM 
    Tn_star(:,tt)  = diag(e_star' * A * e_star / n^2);  % bootstrap statistics  CvM      
    TnU_star(:,tt) = StatCvMRand(hX, theta1, e_star);   %                       CvM: uniform approximation
    TnN_star(:,tt) = StatCvMRand(hX, theta2, e_star);   %                       CvM: gaussian approximation
    KnU_star(:,tt) = StatKSRand(hX, theta1, e_star);    %                       KS: uniform approximation
    KnN_star(:,tt) = StatKSRand(hX, theta2, e_star);    %                       KS: gaussian approximation                 
end % tt in 1:N
%% ---------- Compute rejection rates ---------- %%
power(1,:) = RejectRate(Tn,  Tn_star,  alpha);                       % CvM 
power(2,:) = RejectRate(TnU, TnU_star, alpha);                       % CvM: uniform approximation
power(3,:) = RejectRate(TnN, TnN_star, alpha);                       % CvM: gaussian approximation
power(4,:) = RejectRate(KnU, KnU_star, alpha);                       % KS: uniform approximation
power(5,:) = RejectRate(KnN, KnN_star, alpha);                       % KS: gaussian approximation
EstiPara = [mean(Parabeta)-beta'; sum(bsxfun(@minus,Parabeta,beta').^2)/N; std(Parabeta)];  % [bias; mse; std]
end