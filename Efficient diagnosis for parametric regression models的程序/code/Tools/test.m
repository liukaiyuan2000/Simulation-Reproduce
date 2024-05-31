function [power, EstiPara] = test(n, C, setting, alpha, N, B, M)
%%% Compute KS statistic: random approximation, absolute sign is inside the integral; slower
%%%------Test procedure------%%% 
% input:  n        1*1   sample size
%         C        5*1   distance between HO and H1
%         setting  1*1   select model
%         alpha    1*int significant level 
%         N        1*1   number of trials
%         B        1*1   number of bootstrap: determine critical values 
%         M        1*1   number of resampling: estimate projection 
% output: power    7*int empirical rejection rates
%         EstiPara cell  
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
Zn  = zeros(1,N); Zn_star  = zeros(B,N); Parabeta_Z = zeros(N,p); hX_Z = zeros(n,p); % Zhang(2015)
Jn  = zeros(1,N);                        Parabeta_J = zeros(N,p); hX_J = zeros(n,p); % Zhao(2018)  
%%
for tt = 1 : N
    tt
    %% ---------- Generate data ---------- %%
    [tildeX, tildeY, Z, U] = GeneData(n, C, beta, setting);
    %% ------------ Error correction ------------ %%
    UU = repmat(U,1,n) - (repmat(U,1,n))';
    %------Local linear------%
    h = 2.34 * std(U) * n^(-1/3);  
    Kh = 3/4 * (1-(UU/h).^2) .* (abs(UU)<=h) / h; % standardized Ep kernel function
    hY = CaliLocaLine(tildeY, UU, Kh); 
    for i = 1 : size(tildeX,2)
        hX(:,i) = CaliLocaLine(tildeX(:,i), UU, Kh);
    end
    if size(Z, 2) > 0
        hX(:,p-size(Z,2)+1:p) = Z;
    end
    %---Zhang(2015): local constant---%
    d =  std(U) * n^(-1/3); 
    Kd = 3/4 * (1-(UU/d).^2) .* (abs(UU/d)<=1) / d; 
    hY_Z = CaliLocaCons(tildeY, Kd, 0);
    for i = 1 : size(tildeX,2)
        hX_Z(:,i) = CaliLocaCons(tildeX(:,i), Kd, 0);
    end
    if size(Z, 2) > 0
        hX_Z(:,p-size(Z,2)+1:p) = Z;
    end
    %---Zhao(2018): local constant---%
    g = 0.5 * n^(-1/7); 
    Kg = (15/8) * (3/4) * (1-(7/3)*(UU/g).^2) .* (1-(UU/g).^2) .* (abs(UU/g)<=1) / g;  % standardized 4th kernel function
    hY_J  = CaliLocaCons(tildeY, Kg, 1);
    for i = 1 : size(tildeX,2)
        hX_J(:,i) = CaliLocaCons(tildeX(:,i), Kg, 1);
    end
    if size(Z, 2) > 0
        hX_J(:,p-size(Z,2)+1:p) = Z;
    end
    %% ---- Estimator of beta: (N)LSE ---- %%       
    [Parabeta(tt,:),   e  ] = EstiBeta(hX,   hY,   setting); % this paper
    [Parabeta_Z(tt,:), e_Z] = EstiBeta(hX_Z, hY_Z, setting); % Zhang(2015) 
    [Parabeta_J(tt,:), e_J] = EstiBeta(hX_J, hY_J, setting); % Zhao(2018) 
    %% ---------- Compute test statistics ---------- %%
    theta1 = unifrnd(-1,1,p,M); 
    theta2 = mvnrnd(zeros(M,1),eye(M),p); 
    l_X = exp(0.5*hX_Z(:,2)) - mean(exp(0.5*hX_Z(:,2))); 
    [Tn(tt), A] = StatCvM(hX, e);                       % CvM 
    TnU(tt) = StatCvMRand(hX, theta1, e);               % CvM: uniform approximation
    TnN(tt) = StatCvMRand(hX, theta2, e);               % CvM: gaussian approximation
    KnU(tt) = StatKSRand( hX, theta1, e);               % KS: uniform approximation
    KnN(tt) = StatKSRand( hX, theta2, e);               % KS: gaussian approximation
    Zn(tt)  = StatZhang(l_X, e_Z);                      % Zhang(2015)         
    Jn(tt)  = StatZhao (hX_J, e_J);                     % Zhao(2018)    
    %% --------- Bootstrap resampling --------- %%
    r = unifrnd(0,1,n,B);
    V = (1+sqrt(5))/2 * (r<(5-sqrt(5))/10) + (1-sqrt(5))/2 * (r>=(5-sqrt(5))/10);
    e_star   = BootResi(hX,   hY,   e,   V, setting);   % bootstrap residuals   CvM 
    e_star_Z = BootResi(hX_Z, hY_Z, e_Z, V, setting);   %                       Zhang(2015)
    Tn_star(:,tt)  = diag(e_star' * A * e_star / n^2);  % bootstrap statistics  CvM      
    TnU_star(:,tt) = StatCvMRand(hX, theta1, e_star);   %                       CvM: uniform approximation
    TnN_star(:,tt) = StatCvMRand(hX, theta2, e_star);   %                       CvM: gaussian approximation
    KnU_star(:,tt) = StatKSRand( hX, theta1, e_star);   %                       KS: uniform approximation
    KnN_star(:,tt) = StatKSRand( hX, theta2, e_star);   %                       KS: gaussian approximation          
    for l = 1 : B
        Zn_star(l,tt) = StatZhang(l_X, e_star_Z(:,l));  %                       Zhang(2015)
    end       
end % tt in 1:N
%% ---------- Compute rejection rates ---------- %%
z = norminv(1-alpha,0,1);                                            % upper alpha quantile of N(0,1)
power(1,:) = RejectRate(Tn,  Tn_star,  alpha);                       % CvM 
power(2,:) = RejectRate(TnU, TnU_star, alpha);                       % CvM: uniform approximation
power(3,:) = RejectRate(TnN, TnN_star, alpha);                       % CvM: gaussian approximation
power(4,:) = RejectRate(KnU, KnU_star, alpha);                       % KS: uniform approximation
power(5,:) = RejectRate(KnN, KnN_star, alpha);                       % KS: gaussian approximation
power(6,:) = RejectRate(Zn,  Zn_star,  alpha);                       % Zhang(2015)
power(7,:) = mean(repmat(Jn,length(alpha),1) >= repmat(z',1,N), 2);  % Zhao(2018) 
EstiPara{1} = [mean(Parabeta  )-beta'; sum(bsxfun(@minus,Parabeta,  beta').^2)/N; std(Parabeta)];  % [bias; mse; std]
EstiPara{2} = [mean(Parabeta_Z)-beta'; sum(bsxfun(@minus,Parabeta_Z,beta').^2)/N; std(Parabeta_Z)]; 
EstiPara{3} = [mean(Parabeta_J)-beta'; sum(bsxfun(@minus,Parabeta_J,beta').^2)/N; std(Parabeta_J)];
end