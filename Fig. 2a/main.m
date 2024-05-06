clear;clc;close all;
N =1000; % total users
J=1;
K=100;
L_range=2:2:10;
M_range =[15;32;200]; %The number of antennas
ite = 10000; % All experiments are repeated ite times
mode=1;
K_e=zeros(length(M_range), length(L_range), ite);

for l = 1:length(L_range)
    for m = 1:length(M_range)
        N1=N;
        M=M_range(m);
        L=L_range(l);
        fprintf('L=%d, M=%d\n',L,M);
        for j = 1:ite
            % fprintf('L=%d, M=%d, ite=%d\n',L,M,j);
            K_Est=est_test( N1, K, L, J, M,  mode);
            K_e(m, l, j) = abs(K_Est-K);
        end
    end
end

K_e = K_e / K;
K_e_var = var(K_e, 1, 3);
K_e_mean = mean(K_e, 3);


%% error and var

figure;
errorbar(L_range, K_e_mean(1, :), K_e_var(1, :),'LineWidth',2, 'Color', [0.55,0.80,0.84]); hold on;
errorbar(L_range, K_e_mean(2, :), K_e_var(2, :), 'LineWidth',2, 'Color', [0.90,0.52,0.43]); hold on;
errorbar(L_range, K_e_mean(3, :), K_e_var(3, :), 'LineWidth',2, 'Color', [0.69,0.69,0.69]); hold on;
hp1=plot(L_range, NaN(size(L_range)) ,'-x', 'Color', [0.55,0.80,0.84],'LineWidth',2);
hp2 = plot(L_range, NaN(size(L_range)), '-o', 'Color', [0.90,0.52,0.43],'LineWidth',2);
hp3 = plot(L_range, NaN(size(L_range)), 'Color', [0.69,0.69,0.69],'LineWidth',2);
grid on;
ylim([0.02 0.28]);
legend([hp1, hp2, hp3], {'$M=15$', '$M=32$', '$M=200$'}, 'Interpreter', 'latex');
xlabel({'$L_{\mathrm{I}}$'}, 'Interpreter', 'latex');
ylabel({'$E_K$ '}, 'Interpreter', 'latex');
title({'$K_e$ errorbar with respect to $L_{\mathrm{I}}$'}, 'Interpreter', 'latex');
