clear;clc;
N =1000; % total users
J=1;
K=200;
L_det = 177;
M_range =32; %The number of antennas
ite = 10000; % All experiments are repeated ite times
mode=1;%
Num_step=2000;
time_set= 0.01:0.01:0.25;

CDE_new=zeros(length(time_set),Num_step+1);
FAE_new=zeros(length(time_set),Num_step+1);
cov_time=zeros(length(time_set), 1);


for l = 1:length(time_set)
    max_time=time_set(l);
    fprintf('time=%f\n', max_time);
    N1=N;
    M=M_range;
    %channel
    sigma2s = ones(N1,1);  % large-scale fading
    H = channelGeneration(N1,M,sigma2s);
    for j = 1:ite
        % K_est=est_test(H, N, K, L_est, J, M,  mode);
        [act_set,gamma_cd_new, cov_time1]=test(H, N1, K, L_det, J, M, mode,max_time);
        Thd_cd=0:0.0005:1;
        cov_time(l) = cov_time(l)+cov_time1;
        %%calculate CD&FA of cd
        for s=1:Num_step+1
            thd=Thd_cd(s);
            actset_es_new=find(gamma_cd_new>thd);
            corr_detc=intersect(act_set,actset_es_new);
            fs_dec=setdiff(actset_es_new,corr_detc);
            CDE_new(l,s) = CDE_new(l,s) + length(corr_detc)/K;
            FAE_new(l,s) = FAE_new(l,s) + length(fs_dec)/(N1-K);

        end
    end
end

CDE_new=CDE_new/ite;
FAE_new=FAE_new/ite;
MDE_new=1-CDE_new;
cov_time = cov_time / ite;
%draw
n=length(time_set);
m = Num_step + 1;
sorted_indices_new=zeros(n,m);
value_4=zeros(n,1);
for i=1:n
    [~,sorted_indices_new(i,:)]= sort(abs(MDE_new(i,:)-FAE_new(i,:)));
    index_new=sorted_indices_new(i,1);
    value_4(i)=(FAE_new(i,index_new)+MDE_new(i,index_new))/2;
end

figure;
p4 = semilogy(cov_time,value_4,'o-', LineWidth=2);hold on;
legend(p4, {'CD with $K=200, L=177$'}, 'Interpreter', 'latex');
xlabel('CPU time');
ylabel('Detection error probability');
grid on;
