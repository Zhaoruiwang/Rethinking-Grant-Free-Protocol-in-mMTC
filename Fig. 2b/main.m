clear;clc;close all;
N =1000; % total users
J=1;
K_range=50 : 150 : 950;
L_range=4;
M_range =[15;32;200]; %The number of antennas
ite = 100; % All experiments are repeated ite times
mode=1;
% K_data=zeros(length(M_range),length(L_range),ite);
K_mean=zeros(length(K_range),length(M_range));
K_var=zeros(length(K_range),length(M_range));

for k = 1:length(K_range)
    K=K_range(k);
    for m = 1:length(M_range)
        N1=N;
        M=M_range(m);
        L=L_range;
        for j = 1:ite
            fprintf('L=%d, M=%d, ite=%d\n',L,M,j);
            K_Est=est_test( N1, K, L, J, M,  mode);
            delta=abs(K_Est-K) - K_mean(k,m);
            K_mean(k,m)=K_mean(k,m)+delta/j;
        end
    end
end

%%
figure;
plot(K_range,K_mean(: , 1)./K_range','-x','LineWidth', 2, 'MarkerSize',8,'Color', [0.55,0.80,0.84]); hold on;
plot(K_range,K_mean(: , 2)./K_range','-v','LineWidth', 2, 'MarkerSize',8, 'Color', [0.90,0.52,0.43]); hold on;
plot(K_range,K_mean(: , 3)./K_range','LineWidth', 2, 'MarkerSize',8,'Color', [0.69,0.69,0.69]); hold on;
xlabel('$K$', 'Interpreter', 'latex','FontName','Times New Roman');
ylabel('$E_k$', 'Interpreter', 'latex','FontName','Times New Roman');
grid on;
legend({'$M=15$', '$M=32$', '$M=200$'}, 'Interpreter', 'latex');
%title('Average $K_e$ with respect to $K$, $L1=4$', 'Interpreter', 'latex');
xlim([50 950]);
ylim([0.04 0.22]);
