function K_Est=est_test(N, K, L, J, M,mode)
txPower = 23; % dBm
noisePower = -99; % dBm
txPowerN = 0;
noisePowerN = noisePower + 128 - txPower;
matrx_Type = 'Gaussian';

sigma2s = ones(N,1);  % large-scale fading

% Sequence generation
A = identicalsensingMatrixDesign(N,J,L,matrx_Type,mode);

%Gaussian channel
H = channelGeneration(N,M,sigma2s);

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