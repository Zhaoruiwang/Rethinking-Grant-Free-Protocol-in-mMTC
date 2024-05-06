clear;clc;close all;
N =1000; % total users
J=1;
K=100;
L = 100;
L_est=4;
M_range =15; %The number of antennas
ite = 50; % All experiments are repeated ite times
thd_CD=0.2; % for CD, thd=0.1 for all cases; for benchmark, K=50 thd=0.13; K=90 thd=0.17; K=120 thd=0.2.
mode=1;%
Num_step=200;
CDE_new=zeros(length(M_range),Num_step+1);
FAE_new=zeros(length(M_range),Num_step+1);
CDE_ide=zeros(length(M_range),Num_step+1);
FAE_ide=zeros(length(M_range),Num_step+1);
% CDE_per=zeros(length(M_range),Num_step+1);
% FAE_per=zeros(length(M_range),Num_step+1);
CDE_act=zeros(length(M_range),Num_step+1);
FAE_act=zeros(length(M_range),Num_step+1);
CDE_sel=zeros(length(M_range),Num_step+1);
FAE_sel=zeros(length(M_range),Num_step+1);


cov_time=zeros(4, 1);
iter_time=zeros(4, 1);
Coordinates=zeros(4, 1);
NumberofGradients=zeros(4, 1);

for l = 1:length(M_range)
    N1=N;
    M=M_range(l);
    %channel
    for j = 1:ite
        j
        % fprintf('ite=%d\n',j);
        sigma2s = ones(N1,1);  % large-scale fading
        H = channelGeneration(N1,M,sigma2s);
        %estimate K
        K_est=est_test(H,N, K, L_est, J, M,mode);
        %activity detection
        [cov_time1, iter_time1,Coordinates1,NumberofGradients1,act_set,gamma_cd_new,gamma_cd_ide,gamma_cd_act,gamma_cd_sel] = test(H , N1, K, K_est, L, J, M, thd_CD, mode);
        iter_time=iter_time+iter_time1;
        cov_time=cov_time+cov_time1;
        Coordinates=Coordinates+Coordinates1;
        NumberofGradients=NumberofGradients+NumberofGradients1;
        Thd_cd=[0:0.005:1];
        %%calculate CD&FA of cd
        for s=1:Num_step+1

            thd=Thd_cd(s);

            actset_es_new=find(gamma_cd_new>thd);
            corr_detc=intersect(act_set,actset_es_new);
            fs_dec=setdiff(actset_es_new,corr_detc);
            CDE_new(l,s) = CDE_new(l,s) + length(corr_detc)/K;
            FAE_new(l,s) = FAE_new(l,s) + length(fs_dec)/(N1-K);

            actset_es_ide=find(gamma_cd_ide>thd);
            corr_detc=intersect(act_set,actset_es_ide);
            fs_dec=setdiff(actset_es_ide,corr_detc);
            CDE_ide(l,s) = CDE_ide(l,s) + length(corr_detc)/K;
            FAE_ide(l,s) = FAE_ide(l,s) + length(fs_dec)/(N1-K);

            actset_es_act=find(gamma_cd_act>thd);
            corr_detc=intersect(act_set,actset_es_act);
            fs_dec=setdiff(actset_es_act,corr_detc);
            CDE_act(l,s) = CDE_act(l,s) + length(corr_detc)/K;
            FAE_act(l,s) = FAE_act(l,s) + length(fs_dec)/(N1-K);

            actset_es_sel=find(gamma_cd_sel>thd);
            corr_detc=intersect(act_set,actset_es_sel);
            fs_dec=setdiff(actset_es_sel,corr_detc);
            CDE_sel(l,s) = CDE_sel(l,s) + length(corr_detc)/K;
            FAE_sel(l,s) = FAE_sel(l,s) + length(fs_dec)/(N1-K);

        end
    end
end

CDE_new=CDE_new/ite;
FAE_new=FAE_new/ite;
CDE_ide=CDE_ide/ite;
FAE_ide=FAE_ide/ite;
CDE_act=CDE_act/ite;
FAE_act=FAE_act/ite;
CDE_sel=CDE_sel/ite;
FAE_sel=FAE_sel/ite;

MDE_new=1-CDE_new;
MDE_ide=1-CDE_ide;
MDE_act=1-CDE_act;
MDE_sel = 1-CDE_sel;

NumberofGradients=NumberofGradients/ite;
Coordinates=Coordinates/ite;
iter_time=iter_time/ite;
cov_time=cov_time/ite;

% %dataset = [Coordinates'; NumberofGradients'];
%
% groups=['Random CD', 'Ideal CD', 'active CD''FCD:D=2,a=0.01'];
% data=[Coordinates(1) NumberofGradients(1);Coordinates(2) NumberofGradients(2);Coordinates(3) NumberofGradients(3);Coordinates(4) NumberofGradients(4)];
%
% figure;
% bar(data);
%

NumberofGradients(1)=0;
NumberofGradients(2)=0;
x=[1 2 3 4];
data=[Coordinates(1) NumberofGradients(1);Coordinates(2) NumberofGradients(2);Coordinates(3) NumberofGradients(3);Coordinates(4) NumberofGradients(4)];
data1=round(Coordinates,1);
data2=round(NumberofGradients,1);
X=1:size(data,1);
bar(data,'grouped','BarWidth',1);

for i=1:length(x)
    text(x(i)-0.155,data1(i)+1,num2str(data1(i)),'HorizontalAlignment','center','VerticalAlignment','bottom');
    text(x(i)+0.155,data2(i)+1,num2str(data2(i)),'HorizontalAlignment','center','VerticalAlignment','bottom');
end

algorithmNames = {'CD', 'Ideal CD','Active CD','K-CD'};
metricNames = {'Number of coordinates', 'Number of gradients'};
%bar(x,dataset', 'grouped');
legend(metricNames, 'Location', 'northeast', 'Orientation', 'vertical');
% xlabel('Algorithms');
% ylabel('Number of coordinates/gradients');
set(gca, 'xticklabel', algorithmNames,'YGrid', 'on');
axis([0.5,4.5,0,32000]);
% %axis padded;
%
% ColumnNames={'Random CD', 'Ideal CD', 'Only-a-Few CD'}
% NumofCoorperI=round(Coordinates./iter_time,1);
% NumofGradperI=round(NumberofGradients./iter_time,1);
% data=[ColumnNames',num2cell(round(iter_time,1)),num2cell(NumofCoorperI),num2cell(NumofGradperI)];
% f=figure;
% uitable(f,'Data',data,'ColumnName',{'Algorithms','iteration times','Coordinates per iteration','Gradients per iteration'});