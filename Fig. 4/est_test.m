function K_Est=est_test(H,N, K, L, J, M,mode)

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

txPowerN = 0;
noisePowerN = noisePower + 128 - txPower;
matrx_Type = 'Gaussian';
% large-scale fading
% sigma2s = ones(N,1);  

% Sequence generation
A = identicalsensingMatrixDesign(N,J,L,matrx_Type,mode);

%Gaussian channel
% H = channelGeneration(N,M,sigma2s);

% Sparse signal
[x,user_supp,supp,user_idx,data_idx,sigma2n] = signalGeneration(N,K,L,J,M,H,txPowerN,noisePowerN,mode);
act_set=find(supp);
% Additive noise
w = sqrt(1/2)*(randn(L,M)+1i*randn(L,M))*sqrt(sigma2n);

% System model
y = A * x + w;
sampCov = (1/M)*(y*y'); sigma2 = sigma2n;

K_Est=Estimate_K(A,sigma2,sampCov);
K_Est=round(K_Est);

end