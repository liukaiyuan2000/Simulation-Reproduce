clear
clc
%% ---------- Add path ---------- %%
addpath([pwd, '\Tools\']); % Function files in the specified path can be called
%% ---------- Simulation settings ---------- %%
N = 500;                        % number of trails
B = 1000;                        % number of bootstrap: compute critical values 
M = 100;                        % number of resampling: estimate projection 
alpha = [0.01,0.025,0.05,0.10]; % nominal level
hvalue = (1.34:0.2:3.34);
cvalue = (0.0 : 0.1 : 0.4); % Distance between H0 and H1n
for setting = [1 3 5]      % select model
    setting
    for n = [100 300]     % sample size
        n
        %% -------- Calculate and save the results -------- %%  
        filename = ['DebugHn_N',int2str(N), '_B',int2str(B), '_setting',int2str(setting), '_n',int2str(n)];
        fid = fopen([filename,'.txt'], 'w');
        
        EstiPara = cell(length(hvalue),length(cvalue));
        for ii = 1 : length(hvalue)
            hn = hvalue(ii)
            for c = 1 : length(cvalue)
                C = cvalue(c)
                tic
                [power, EstiPara{ii,c}] = DebugHnCvM(n, C, setting, alpha, N, B, hn);
                %[power, EstiPara{ii,c}] = DebugHn(n, C, setting, alpha, N, B, M, hn);
                toc
                fprintf(fid,'%s %d %s %.2f %s %.2f %s\t','------ n=',n,'C=',cvalue(c),'hn=',hvalue(ii),'------');
                fprintf(fid,'\n');
                for i = 1 : size(power,1)
                    fprintf(fid,'%.4f\t',power(i,:)); % output: according to the line
                    fprintf(fid,'\n');
                end
                fprintf(fid,'\n');
            end % C
        end % hn
        fclose(fid);
%         save([filename,'_EstiPara.mat'], 'EstiPara')
    end
end
