clear
clc
%% ---------- Add path ---------- %%
addpath([pwd, '\Tools\']); % Function files in the specified path can be called
%% ---------- Simulation settings ---------- %%
N = 500;                        % number of trails
B = 1000;                       % number of bootstrap: compute critical values 
M = 50;                         % number of resampling: estimate projection 
alpha = [0.01,0.025,0.05,0.10]; % nominal level
cvalue = (0 : 0.1 : 0.4);%      % Distance between H0 and H1n 
setting = [1,2,3,4,5];   % select model
n = [100]; % sample size
[I_setting,J_n] = meshgrid(setting,n);
for i = 1 : numel(I_setting)
    %% -------- Calculate and save the results -------- %%  
    filename = ['1output_N',int2str(N), '_B',int2str(B), '_setting',int2str(I_setting(i)), '_n',int2str(J_n(i))];
    fid = fopen([filename,'.txt'], 'w');
    for c = 1 : length(cvalue)
        C = cvalue(c)
        tic
        [power, ~] = test(J_n(i), C, I_setting(i), alpha, N, B, M); 
        toc
        fprintf(fid,'%s %d %s %.2f %s\t','------ n=',J_n(i),'C=',cvalue(c),'------');
        fprintf(fid,'\n');
        for j = 1 : size(power,1)
            fprintf(fid,'%.4f\t',power(j,:)); % output: according to the line
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end % C
    fclose(fid);    
end
