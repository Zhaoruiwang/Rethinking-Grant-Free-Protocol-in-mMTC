clear;clc;close all;
N =1000; % total users
J=1;
K=100;
L_est = 4;
M_range =32; %The number of antennas
ite = 10000; % All experiments are repeated ite times
mode=1;%
Num_step=2000;
time_set= 0.01:0.01:0.25;
L_sum =  zeros(length(time_set), 1);
CDE_sel=zeros(length(time_set),Num_step+1);
FAE_sel=zeros(length(time_set),Num_step+1);
cov_time_100 = zeros(length(time_set), 1);

for l = 1:length(time_set)
    max_time=time_set(l);
    fprintf('time=%f\n',time_set(l));
    N1=N;
    M=M_range;
    %channel
    sigma2s = ones(N1,1);  % large-scale fading
    H = channelGeneration(N1,M,sigma2s);
    for j = 1:ite
        K_est=est_test(H, N, K, L_est, J, M,  mode);
        L_det = L_allocation(K_est);
        L_sum(l) = L_sum(l) + L_det;
        [act_set,gamma_cd_sel, cov_time1]=test(H, N1, K, K_est, L_det, J, M, mode,max_time);
        cov_time_100(l) = cov_time_100(l) + cov_time1;
        Thd_cd=0:0.0005:1;
        %%calculate CD&FA of cd
        for s=1:Num_step+1
            thd=Thd_cd(s);
            actset_es_sel=find(gamma_cd_sel>thd);
            corr_detc=intersect(act_set,actset_es_sel);
            fs_dec=setdiff(actset_es_sel,corr_detc);
            CDE_sel(l,s) = CDE_sel(l,s) + length(corr_detc)/K;
            FAE_sel(l,s) = FAE_sel(l,s) + length(fs_dec)/(N1-K);

        end
    end
end

CDE_sel=CDE_sel/ite;
FAE_sel=FAE_sel/ite;
MDE_sel=1-CDE_sel;
cov_time_100 = cov_time_100 / ite;
L_sum = L_sum / ite;
n=length(time_set);
m = Num_step + 1;
sorted_indices_sel=zeros(n,m);
value_100=zeros(n,1);
for i=1:n
    [~,sorted_indices_sel(i,:)]= sort(abs(MDE_sel(i,:)-FAE_sel(i,:)));
    index_sel=sorted_indices_sel(i,1);
    value_100(i)=(FAE_sel(i,index_sel)+MDE_sel(i,index_sel))/2;
end

figure;
p1 = semilogy(cov_time_100,value_100,'o-', LineWidth=2);hold on;
legend(p1, {'$K=100,L_{I}=4, L_{II}=83, D=2,\alpha=0.01$'}, 'Interpreter', 'latex');
xlabel('CPU time');
ylabel('Detection error probability');
grid on;
