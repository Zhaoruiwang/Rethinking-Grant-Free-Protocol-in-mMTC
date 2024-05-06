function [x,user_supp,supp,user_idx,data_idx,sigma2n]=identicalsignalGeneration(N,K,L,J,M,H,txPowerMax,noisePower,mode)
% activity pattern
user_idx = randperm(N);
user_idx(1:K) = sort(user_idx(1:K)); % user_supp(1:K), the index of active users
user_supp = zeros(N,1);
user_supp(user_idx(1:K)) = 1;

%signal generation
%mode 1: only have one sequence for one user
if mode==1
    data_idx=0;
    supp=zeros(N,1);
    supp(user_idx(1:K))=1;
    x=zeros(N,M);
    x((supp==1),:)= H((supp==1),:);
end

%mode 2: one user have 2^J possible sequence
% Data of active users 
if mode==2
    data_idx = randi(2^J,K,1); % the data indices for active users
    % Combined support
    supp = zeros(2^J*N,1);   
    supp(sort(user_idx(1:K)-1)'*2^J + data_idx) = 1;
    % Signal generation 
    Ne = N*2^J;
    Heff = repelem(H,2^J,1); % eff. channel; channel for sequences of one users are equal.
    x = zeros(Ne,M);
    x((supp==1),:) = Heff((supp==1),:); 
end


%noise setup with power control
sigma2n = (10^((noisePower)/10));
txPower = 10^(txPowerMax/10);
sigma2n = sigma2n/txPower;

end



