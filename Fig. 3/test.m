function [cov_time, act_set, gamma_cd_new, gamma_cd_idc, gamma_cd_act, gamma_cd_sel1,gamma_cd_sel2,gamma_cd_sel3]=test(H,N, K, K_est, L, L_det, J, M,mode,max_time)


txPower = 23; % dBm
noisePower = -99; % dBm

sigma2s = ones(N,1);  % large-scale fading
txPowerN = 0;
noisePowerN = noisePower + 128 - txPower;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrx_Type = 'Gaussian';

% Sequence generation
A = sensingMatrixDesign(N,J,L,matrx_Type,mode);
% A = A(1 : L_det, :);
% Gaussian channel
%H = channelGeneration(N,M,sigma2s);

% Sparse signal
[x,user_supp,supp,user_idx,data_idx,sigma2n] = signalGeneration(N,K,L,J,M,H,txPowerN,noisePowerN,mode);
act_set=find(supp);
% Additive noise
w = sqrt(1/2)*(randn(L,M)+1i*randn(L,M))*sqrt(sigma2n);
% System model
y = A * x + w;

sampCov = (1/M)*(y*y'); sigma2 = sigma2n;


gamma = zeros(N,1);

[gamma_cd_new,  ~, cov_time_new, ~, ~] = Random_CD_new_time(gamma, A, sampCov, sigma2, J,N, max_time);

[gamma_cd_act,~, cov_time_act, ~] = active_CD_time(gamma, A, sampCov, sigma2, J,N, max_time);

[~,  act_set_es_new, ~, ~, ~] = Random_CD_new_time_supp(gamma, A, sampCov, sigma2, J,N);

[gamma_cd_idc,~,cov_time_idc, ~, ~]=Ideal_CD_cdbased_time(gamma,A, sampCov, sigma2, act_set_es_new, N, max_time);

[gamma_cd_sel1, ~, cov_time_sel1, ~, ~] = CD_select2_time(gamma, A, sampCov, sigma2, J, N, K_est, 0.1, 0, max_time);

[gamma_cd_sel2, ~, cov_time_sel2, ~, ~]= CD_select2_time(gamma, A, sampCov, sigma2, J, N, K_est, 0.01, 0, max_time);

[gamma_cd_sel3, ~, cov_time_sel3, ~, ~] = CD_select2_time(gamma, A, sampCov, sigma2, J, N, K_est, 0.01, 2, max_time);

cov_time =  [cov_time_new, cov_time_idc, cov_time_act, cov_time_sel1, cov_time_sel2, cov_time_sel3];


end