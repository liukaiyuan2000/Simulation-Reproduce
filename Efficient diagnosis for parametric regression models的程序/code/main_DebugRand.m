clear
clc
%% ---------- 添加路径 ---------- %%
addpath([pwd, '\Tools\']); % 可调用指定路径中的函数文件
%% ---------- 模拟设置 ---------- %%
N = 500; % 模拟次数
B = 1000; % bootstrap次数 estimate critical value
% M = [20 25 40 50 60 75 80 100 150 200 250 300 350 400 450 500]; 
M = [25 50 75 100 125 150 175 200 225 250 275 300]; % bootstrap次数 estimate projection parameter
alpha = [0.01,0.025,0.05,0.10]; % 显著性水平
n = 100;
for setting = [1 5]
    setting
    filenameU_CvM = ['DebugRand_Uniform_CvM_N',int2str(N), '_B',int2str(B), '_setting',int2str(setting), '_n',int2str(n)];
    filenameN_CvM = ['DebugRand_Normal_CvM_N',int2str(N), '_B',int2str(B), '_setting',int2str(setting), '_n',int2str(n)];
    filenameU_KS = ['DebugRand_Uniform_KS_N',int2str(N), '_B',int2str(B), '_setting',int2str(setting), '_n',int2str(n)];
    filenameN_KS = ['DebugRand_Normal_KS_N',int2str(N), '_B',int2str(B), '_setting',int2str(setting), '_n',int2str(n)];
    fidU_CvM = fopen([filenameU_CvM,'.txt'], 'w');
    fidN_CvM = fopen([filenameN_CvM,'.txt'], 'w');
    fidU_KS = fopen([filenameU_KS,'.txt'], 'w');
    fidN_KS = fopen([filenameN_KS,'.txt'], 'w');
    %% ---------- 变量初始化 ---------- %%
    switch setting
        case 1   % Zhang's example 1-model 1: LM p=2 q=0 
            beta = [2;3];
            cvalue = (0 : 0.1 : 0.4); %(0 : 0.2 : 0.8); 
        case 2   % Zhang's example 1-model 2: NLM p=2 q=0 
            beta = [1;2]; 
            cvalue = (0 : 0.1 : 0.4); 
        case 3   % LM p=5 q=0
            beta = [1;1;1;1;1];
            cvalue = (0 : 0.1 : 0.4); 
        case 4   % NLM p=3 q=2  
            beta = [1;1;-1;-1;-1];
            cvalue = (0 : 0.1 : 0.4);
        case 5   % LM p=6 q=4
            beta = [1;1;1;1;-1;-1;-1;-1;-1;-1];
            cvalue = (0 : 0.1 : 0.4);
    end 
    p = length(beta); 
    tXtemp = cell(N,1);
    tYtemp = cell(N,1);
    Ztemp = cell(N,1);
    Utemp = cell(N,1);
    hX = zeros(n,p); 
    %--检验统计量-----bootstrap检验统计量---------参数估计值----------%
    Tn = zeros(1,N); Tn_star = zeros(B,N); Parabeta   = zeros(N,p); % 本文CvM 
    Tn_estiU = cell(length(M),1);
    Tn_esti_starU = cell(length(M),1);
    Kn_estiU = cell(length(M),1);
    Kn_esti_starU = cell(length(M),1);
    Tn_estiN = cell(length(M),1);
    Tn_esti_starN = cell(length(M),1);
    Kn_estiN = cell(length(M),1);
    Kn_esti_starN = cell(length(M),1);
    for i = 1:length(M)
        Tn_estiU{i} = zeros(1,N);
        Tn_esti_starU{i} = zeros(B,N);
        Kn_estiU{i} = zeros(1,N); 
        Kn_esti_starU{i} = zeros(B,N);
        Tn_estiN{i} = zeros(1,N);
        Tn_esti_starN{i} = zeros(B,N);
        Kn_estiN{i} = zeros(1,N); 
        Kn_esti_starN{i} = zeros(B,N);
    end

    for c = 1 : length(cvalue)
        tic
        C = cvalue(c)
        for tt = 1 : N
        %     tt
            %% ---------- 产生模拟数据 ---------- %%
             [tildeX, tildeY, Z, U] = GeneData(n, C, beta, setting);
             tXtemp{tt} = tildeX;
             tYtemp{tt} = tildeY;
             Ztemp{tt} = Z;
             Utemp{tt} = U;
            %% ---------- 误差校正 ---------- %%
            UU = repmat(U,1,n) - (repmat(U,1,n))';
            %---本文：局部线性---%
            h = 2.34 * std(U) * n^(-1/3);  
            Kh = 3/4 * (1-(UU/h).^2) .* (abs(UU)<=h) / h; % standardized Ep kernel function
            hY = CaliLocaLine(tildeY, UU, Kh); 
            for i = 1 : size(tildeX,2)
                hX(:,i) = CaliLocaLine(tildeX(:,i), UU, Kh);
            end
            if size(Z, 2) > 0
                hX(:,p-size(Z,2)+1:p) = Z;
            end
            %% ---------- 参数估计:（非）线性最小二乘 ---------- %%       
            [Parabeta(tt,:), e] = EstiBeta(hX, hY, setting); % 本文
            %% ---------- 计算检验统计量 ---------- %%
            thetaU = cell(length(M),1);
            thetaN = cell(length(M),1);
            for i = 1:length(M)
                thetaU{i} = unifrnd(-1,1,p,M(i));      
                thetaN{i} = mvnrnd(zeros(M(i),1),eye(M(i)),p); 
                Tn_estiU{i}(tt) = StatCvMRand(hX, thetaU{i}, e);         
                Kn_estiU{i}(tt) = StatKSRand1(hX, thetaU{i}, e);
                Tn_estiN{i}(tt) = StatCvMRand(hX, thetaN{i}, e);         
                Kn_estiN{i}(tt) = StatKSRand1(hX, thetaN{i}, e);              
            end    
            [Tn(tt), A] = StatCvM(hX, e);               
            %% ---------- bootstrap重抽样 ---------- %%
            r = unifrnd(0,1,n,B);
            V = (1+sqrt(5))/2 * (r<(5-sqrt(5))/10) + (1-sqrt(5))/2 * (r>=(5-sqrt(5))/10);
            e_star = BootResi(hX, hY, e, V, setting); % bootstrap残差
            %---bootstrap统计量---%
            Tn_star(:,tt) = diag(e_star' * A * e_star / n^2);  
            for i = 1:length(M)
                Tn_esti_starU{i}(:,tt) = StatCvMRand(hX, thetaU{i}, e_star);    
                Kn_esti_starU{i}(:,tt) = StatKSRand1(hX, thetaU{i}, e_star); 
                Tn_esti_starN{i}(:,tt) = StatCvMRand(hX, thetaN{i}, e_star);    
                Kn_esti_starN{i}(:,tt) = StatKSRand1(hX, thetaN{i}, e_star);             
            end 
        end % tt in 1:N
        %% ---------- 计算拒绝率 ---------- %%
        para = [mean(Parabeta)-beta'; sum(bsxfun(@minus,Parabeta,beta').^2)/N; std(Parabeta)]
        powerU_CvM = zeros(length(M)+1,length(alpha));
        powerN_CvM = zeros(length(M)+1,length(alpha));
        powerU_KS = zeros(length(M)+1,length(alpha));
        powerN_KS = zeros(length(M)+1,length(alpha));

        powerU_CvM(1,:) = RejectRate(Tn, Tn_star, alpha);
        powerN_CvM(1,:) = powerU_CvM(1,:);
        powerU_KS(1,:) = powerU_CvM(1,:);
        powerN_KS(1,:) = powerU_CvM(1,:);

        for i = 1 : length(M)
            powerU_CvM(i+1,:) = RejectRate(Tn_estiU{i}, Tn_esti_starU{i}, alpha);
            powerN_CvM(i+1,:) = RejectRate(Tn_estiN{i}, Tn_esti_starN{i}, alpha);
            powerU_KS(i+1,:) = RejectRate(Kn_estiU{i}, Kn_esti_starU{i}, alpha);
            powerN_KS(i+1,:) = RejectRate(Kn_estiN{i}, Kn_esti_starN{i}, alpha);        
        end
        %% 输出结果   
        fprintf(fidU_CvM,'%s %d %s %.2f %s\t','------ n=',n,'C=',cvalue(c),'------');
        fprintf(fidU_CvM,'\n');
        for i = 1 : length(M)+1
            fprintf(fidU_CvM,'%.4f\t',powerU_CvM(i,:)); % output: according to the line
            fprintf(fidU_CvM,'\n');
        end
        fprintf(fidU_CvM,'\n');

        fprintf(fidN_CvM,'%s %d %s %.2f %s\t','------ n=',n,'C=',cvalue(c),'------');
        fprintf(fidN_CvM,'\n');
        for i = 1 : length(M)+1
            fprintf(fidN_CvM,'%.4f\t',powerN_CvM(i,:)); % output: according to the line
            fprintf(fidN_CvM,'\n');
        end
        fprintf(fidN_CvM,'\n');

        fprintf(fidU_KS,'%s %d %s %.2f %s\t','------ n=',n,'C=',cvalue(c),'------');
        fprintf(fidU_KS,'\n');
        for i = 1 : length(M)
            fprintf(fidU_KS,'%.4f\t',powerU_KS(i,:)); % output: according to the line
            fprintf(fidU_KS,'\n');
        end
        fprintf(fidU_KS,'\n');

        fprintf(fidN_KS,'%s %d %s %.2f %s\t','------ n=',n,'C=',cvalue(c),'------');
        fprintf(fidN_KS,'\n');
        for i = 1 : length(M)
            fprintf(fidN_KS,'%.4f\t',powerN_KS(i,:)); % output: according to the line
            fprintf(fidN_KS,'\n');
        end
        fprintf(fidN_KS,'\n');
        toc
    end % cvalue
    
    fclose(fidU_CvM);  
    fclose(fidN_CvM);  
    fclose(fidU_KS);  
    fclose(fidN_KS);  
end % setting
