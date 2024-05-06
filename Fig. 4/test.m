function [cov_time,iter_time,Coordinates,NumberofGradients,act_set,gamma_cd_new,gamma_cd_ide,gamma_cd_act,gamma_cd_sel]=test(H,N, K, K_est, L, J, M,thd_CD,mode)


txPower = 23; % dBm
noisePower = -99; % dBm

sigma2s = ones(N,1);  % large-scale fading
txPowerN = 0;
noisePowerN = noisePower + 128 - txPower;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrx_Type = 'Gaussian';

% Sequence generation
A = sensingMatrixDesign(N,J,L,matrx_Type,mode); %L*N
%A_truncated = A(1 : L_det, :);
% Gaussian channel
%H = channelGeneration(N,M,sigma2s);

% Sparse signal
[x,user_supp,supp,user_idx,data_idx,sigma2n] = signalGeneration(N,K,L,J,M,H,txPowerN,noisePowerN,mode); %x is N * M
act_set=find(supp);
% Additive noise
w = sqrt(1/2)*(randn(L,M)+1i*randn(L,M))*sqrt(sigma2n);
%w_truncated = sqrt(1/2)*(randn(L_det, M)+1i*randn(L_det, M))*sqrt(sigma2n);
% System model
y = A * x + w; % L * M
%y_truncated = A_truncated * x + w_truncated;
sampCov = (1/M)*(y*y'); sigma2 = sigma2n;
%sampCov_truncated = (1 / M) * (y_truncated * y_truncated');
gamma = zeros(N,1);

[gamma_cd_new,  act_set_es_new,cov_time_new,iter_time_new,Coord_total_new,NumberofGradients_new] = Random_CD_new(gamma,A, sampCov, sigma2, J,thd_CD, N);

[gamma_cd_ide, ~, cov_time_ide,iter_time_ide, Coord_total_ide, NumberofGradients_ide] = Ideal_CD_cdbased(gamma, A, sampCov, sigma2, act_set_es_new, N);

[gamma_cd_act,~, cov_time_act, iter_time_act, Coord_total_act, NumberofGradients_act] = active_CD(gamma, A, sampCov, sigma2, J,N);

[gamma_cd_sel, ~, cov_time_sel, iter_time_sel, Coord_total_sel, NumberofGradients_sel] = CD_select2(gamma, A, sampCov, sigma2, J, thd_CD, N, K_est, 0.01, 2);

NumberofGradients=[NumberofGradients_new; NumberofGradients_ide; NumberofGradients_act; NumberofGradients_sel];
Coordinates=[Coord_total_new; Coord_total_ide;  Coord_total_act; Coord_total_sel];
iter_time=[iter_time_new; iter_time_ide; iter_time_act; iter_time_sel];
cov_time=[cov_time_new; cov_time_ide; cov_time_act; cov_time_sel];

end