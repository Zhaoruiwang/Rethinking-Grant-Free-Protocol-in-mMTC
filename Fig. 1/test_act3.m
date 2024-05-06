function [ act_set,act_set_es, gamma_cd] = test_act3(H,N, K,  L, J, M,thd_CD,mode,N_iter)
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
epsilon = 1e-3;


% fprintf('CD\n');
[gamma_cd,  act_set_es,cov_time,iter_time,d_cd] = Random_CD_Niter(gamma,A, sampCov, sigma2, J,thd_CD,N,N_iter);
% [gamma_bcd,  act_set_es,cov_time] = Random_BCD(A, sampCov, sigma2, J,thd_BCD);
% [gamma_bcd,  act_set_es,cov_time] = Random_BCD_wp(A, sampCov, sigma2, Q_max,thd_BCD,K);




end