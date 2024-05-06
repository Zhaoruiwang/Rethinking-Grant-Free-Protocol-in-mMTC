function [act_set, gamma_cd_sel, cov_time_sel]=test(H,N, K, K_est, L_det, J, M,mode,max_time)

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
A = sensingMatrixDesign(N,J,L_det,matrx_Type,mode);
% A_truncated = A(1 : L_det, :);
% Gaussian channel
%H = channelGeneration(N,M,sigma2s);

% Sparse signal
[x,user_supp,supp,user_idx,data_idx,sigma2n] = signalGeneration(N,K,L_det,J,M,H,txPowerN,noisePowerN,mode);
act_set=find(supp);
% Additive noise
w = sqrt(1/2)*(randn(L_det, M)+1i*randn(L_det, M))*sqrt(sigma2n);
% w_truncated = sqrt(1/2)*(randn(L_det, M)+1i*randn(L_det, M))*sqrt(sigma2n);
% System model
y = A * x + w;
sampCov = (1/M)*(y*y'); sigma2 = sigma2n;
gamma = zeros(N,1);



[gamma_cd_sel, ~, cov_time_sel, ~, ~] = CD_select2_time(gamma, A, sampCov, sigma2, J, N, K_est, 0.01, 2, max_time);

end