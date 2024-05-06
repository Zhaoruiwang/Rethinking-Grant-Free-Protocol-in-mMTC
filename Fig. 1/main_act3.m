clear;clc;close all;
N = 1000; % total users
J=1;
K_set=[100 200 200];
L_set=[100 100 180];
M_range = 15; %The number of antennas
ite = 10000; % All experiments are repeated ite times
N_iter=200;
thd_CD=0.2; % for CD, thd=0.1 for all cases; for benchmark, K=50 thd=0.13; K=90 thd=0.17; K=120 thd=0.2.
mode=1;%
Num_step=200;
CDE=zeros(length(K_set),Num_step+1);
FAE=zeros(length(K_set),Num_step+1);

N1=N;
M=M_range;
%channel
sigma2s = ones(N1,1);  % large-scale fading
H = channelGeneration(N1,M,sigma2s);
for j = 1:ite
    for l=1:length(K_set)
        K=K_set(l);
        L=L_set(l);
        [ act_set,act_set_es, gamma_cd] = test_act3(H , N1, K, L, J, M, thd_CD, mode,N_iter);
        %%calculate CD&FA of cd
        for s=1:Num_step+1
            Thd_cd=[min(gamma_cd):(max(gamma_cd)-min(gamma_cd))/Num_step:max(gamma_cd)];
            thd=Thd_cd(s);
            actset_es=find(gamma_cd>thd);
            corr_detc=intersect(act_set,actset_es);
            fs_dec=setdiff(actset_es,corr_detc);
            CDE(l,s) = CDE(l,s) + length(corr_detc)/K;
            FAE(l,s) = FAE(l,s) + length(fs_dec)/(N1-K);
        end
    end    
end

CDE=CDE/ite;
FAE=FAE/ite;
MDE=1-CDE;

a=[1:5:200];
b=[1:5:200];
figure;
p1=loglog(MDE(1,a),FAE(1,a),'-o','LineWidth',2); hold on;
p1.Color = "#8BCCD6";
p2=loglog(MDE(2,b),FAE(2,b),'-^','LineWidth',2); hold on;
p2.Color = "#8BCCD6";
p3=loglog(MDE(3,a),FAE(3,a),'-^','LineWidth',2); hold on;
p3.Color = "#E6846E";
p3.LineStyle='- -';
grid on;
xlabel('MDP','FontSize',13);
ylabel('FAP','FontSize',13);
legend({'$K=100, L_{\rm I}=100$','$K=200, L_{\rm I}=100$','$K=200, L_{\rm I}=180$'},'Interpreter','latex','FontSize',14);
axis([1e-5, 1, 1e-5, 1e-1]);