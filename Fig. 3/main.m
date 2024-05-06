clear;clc;
N =1000; % total users
J=1;
K=100;
L=100;
L_est = 4;
L_det = L;
M_range =15; %The number of antennas
ite =  10000; % All experiments are repeated ite times
mode=1;%
Num_step=2000;
time_set= 0.01: 0.01 : 0.1;

CDE_new=zeros(length(time_set),Num_step+1);
FAE_new=zeros(length(time_set),Num_step+1);
CDE_idc=zeros(length(time_set),Num_step+1);
FAE_idc=zeros(length(time_set),Num_step+1);
CDE_sel1=zeros(length(time_set),Num_step+1);
FAE_sel1=zeros(length(time_set),Num_step+1);
CDE_sel2=zeros(length(time_set),Num_step+1);
FAE_sel2=zeros(length(time_set),Num_step+1);
CDE_sel3=zeros(length(time_set),Num_step+1);
FAE_sel3=zeros(length(time_set),Num_step+1);
CDE_sel4=zeros(length(time_set),Num_step+1);
FAE_sel4=zeros(length(time_set),Num_step+1);
CDE_act=zeros(length(time_set),Num_step+1);
FAE_act=zeros(length(time_set),Num_step+1);


cov_time=zeros(length(time_set),6);
iter_time=zeros(length(time_set),6);
Coordinates=zeros(length(time_set),6);

for l = 1:length(time_set)
    max_time=time_set(l);
    fprintf('time=%f\n', max_time);
    N1=N;
    M=M_range;
    %channel
    sigma2s = ones(N1,1);  % large-scale fading
    H = channelGeneration(N1,M,sigma2s);
    for j = 1:ite
        % fprintf('time=%f, ite=%d\n',time_set(l),j);
        K_est=est_test(H, N, K, L_est, J, M,  mode);
        [cov_time1, act_set, gamma_cd_new, gamma_cd_idc, gamma_cd_act ,gamma_cd_sel1,gamma_cd_sel2, gamma_cd_sel3]=test(H , N1, K, K_est, L, L_det, J, M, mode,max_time);
        cov_time(l, :)=cov_time(l, :)+cov_time1;
        Thd_cd=0:0.0005:1;
        %%calculate CD&FA of cd
        for s=1:Num_step+1
            thd=Thd_cd(s);
            actset_es_new=find(gamma_cd_new>thd);
            corr_detc=intersect(act_set,actset_es_new);
            fs_dec=setdiff(actset_es_new,corr_detc);
            CDE_new(l,s) = CDE_new(l,s) + length(corr_detc)/K;
            FAE_new(l,s) = FAE_new(l,s) + length(fs_dec)/(N1-K);

            actset_es_idc=find(gamma_cd_idc>thd);
            corr_detc=intersect(act_set,actset_es_idc);
            fs_dec=setdiff(actset_es_idc,corr_detc);
            CDE_idc(l,s) = CDE_idc(l,s) + length(corr_detc)/K;
            FAE_idc(l,s) = FAE_idc(l,s) + length(fs_dec)/(N1-K);

            actset_es_sel1=find(gamma_cd_sel1>thd);
            corr_detc=intersect(act_set,actset_es_sel1);
            fs_dec=setdiff(actset_es_sel1,corr_detc);
            CDE_sel1(l,s) = CDE_sel1(l,s) + length(corr_detc)/K;
            FAE_sel1(l,s) = FAE_sel1(l,s) + length(fs_dec)/(N1-K);

            actset_es_sel2=find(gamma_cd_sel2>thd);
            corr_detc=intersect(act_set,actset_es_sel2);
            fs_dec=setdiff(actset_es_sel2,corr_detc);
            CDE_sel2(l,s) = CDE_sel2(l,s) + length(corr_detc)/K;
            FAE_sel2(l,s) = FAE_sel2(l,s) + length(fs_dec)/(N1-K);

            actset_es_sel3=find(gamma_cd_sel3>thd);
            corr_detc=intersect(act_set,actset_es_sel3);
            fs_dec=setdiff(actset_es_sel3,corr_detc);
            CDE_sel3(l,s) = CDE_sel3(l,s) + length(corr_detc)/K;
            FAE_sel3(l,s) = FAE_sel3(l,s) + length(fs_dec)/(N1-K);

            actset_es_act=find(gamma_cd_act>thd);
            corr_detc=intersect(act_set,actset_es_act);
            fs_dec=setdiff(actset_es_act,corr_detc);
            CDE_act(l,s) = CDE_act(l,s) + length(corr_detc)/K;
            FAE_act(l,s) = FAE_act(l,s) + length(fs_dec)/(N1-K);
        end
    end
end

cov_time = cov_time / ite;
iter_time = iter_time / ite;

CDE_new=CDE_new/ite;
FAE_new=FAE_new/ite;
CDE_act = CDE_act/ite;
FAE_act = FAE_act /ite;
CDE_idc=CDE_idc/ite;
FAE_idc=FAE_idc/ite;
CDE_sel1=CDE_sel1/ite;
FAE_sel1=FAE_sel1/ite;
CDE_sel2=CDE_sel2/ite;
FAE_sel2=FAE_sel2/ite;
CDE_sel3=CDE_sel3/ite;
FAE_sel3=FAE_sel3/ite;


MDE_new=1-CDE_new;
MDE_idc=1-CDE_idc;
MDE_act = 1 - CDE_act;
MDE_sel1=1-CDE_sel1;
MDE_sel2=1-CDE_sel2;
MDE_sel3=1-CDE_sel3;
