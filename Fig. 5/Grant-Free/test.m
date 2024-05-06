function [act_set, gamma_cd_new, cov_time_CD]=test(H,N, K, L, J, M,mode,max_time)

%,gamma_cd_ideal,gamma_cd_idc] = test(H,N, K, K_est, L, J, M,thd_CD,mode)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this function we perform an experiment on all algorithms
% Note that the power is normalized by 128dBm
% N = 1000;
% K = 150;
% L = 64;
% J = 2;
% M = 128;
% maxN_itera = 50;
txPower = 23; % dBm
noisePower = -99; % dBm

sigma2s = ones(N,1);  % large-scale fading
txPowerN = 0;
noisePowerN = noisePower + 128 - txPower;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrx_Type = 'Gaussian';

% Sequence generation
A = sensingMatrixDesign(N,J,L,matrx_Type,mode);
% A_truncated = A(1 : L_det, :);
% Gaussian channel
%H = channelGeneration(N,M,sigma2s);

% Sparse signal
[x,user_supp,supp,user_idx,data_idx,sigma2n] = signalGeneration(N,K,L,J,M,H,txPowerN,noisePowerN,mode);
act_set=find(supp);
% Additive noise
w = sqrt(1/2)*(randn(L, M)+1i*randn(L, M))*sqrt(sigma2n);
% w_truncated = sqrt(1/2)*(randn(L_det, M)+1i*randn(L_det, M))*sqrt(sigma2n);
% System model
y = A * x + w;
% y_truncated = A_truncated * x + w_truncated;
sampCov = (1/M)*(y*y'); sigma2 = sigma2n;
% sampCov_truncated = (1 / M) * (y_truncated * y_truncated');
gamma = zeros(N,1);

[gamma_cd_new,  ~,cov_time_CD,~,~] = Random_CD_new_time(gamma,A, sampCov, sigma2, J,N,max_time);



end