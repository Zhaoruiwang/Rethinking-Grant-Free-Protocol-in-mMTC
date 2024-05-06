function [y,t,iter_time,NumberofCoord,NumberofGradients] = act_algorithm(A, sampCov, sigma2, flag)

%{
Active set algorithm
flag == 
1: use the original gradient
2: use the scaled gradient (scaled by ||s_n||^2)
%}
 
tic;
[L, N] = size(A);
gamma = zeros(N,1);

epsilon = 1e-3;

NumberofCoord=0;
NumberofGradients=0;

if flag == 1
    scale = qones(N, 1);
elseif flag == 2
    Anorm = sum( A .* conj(A) )';
    scale = Anorm / L;
end

k = 0;

while k < 1000
    
    k = k + 1;
    
    g = grad(gamma, A, sampCov, sigma2);
    NumberofGradients=NumberofGradients+N;
    V_gamma=V(N,gamma,g);

    if norm(V_gamma,'inf') < epsilon
        break
    end
    
    % select the active set
    actset = selectAct(gamma, g, k, scale);
    
    % solve the subproblem on the active set
    [gamma(actset),iter1] = updateAct_SPG2(gamma(actset), A(:, actset), sampCov, sigma2, k, epsilon);
    NumberofCoord=NumberofCoord+length(actset);
    NumberofGradients=NumberofGradients+length(actset)*iter1;
end

iter_time=k;
% fprintf('A total of %d active sets were selected\n', k-1);
y = gamma;
t = toc;
end


function actset = selectAct(gamma, g, k_iter, scale)
%{
This function is used to select a active set
%}
N = length(gamma);
g_scale = g ./ scale .* (gamma == 0);%return 1 if gamma==0
% g1 represents the gradient of the nonpositive coordinates
g_min = min(g_scale);

% c = 0.5;
c = max(0.4, 0.7 - k_iter/20);
% epsilon1 = min(1e-4, 1/10^k_iter);
epsilon1 = 1e-6 * 10^(-k_iter);


% supp = find((gamma == 0 & g < c * g_min & g < 0) | (gamma > 0));
if min(g) > -1
    epsilon2 = 1/10^(k_iter-1);
else
    epsilon2 = inf;
end

while sum(g < -epsilon2) >= N*0.1
    epsilon2 = epsilon2 * 10;
end

while sum(g_scale < c * g_min) >= N*0.1 && c < 1
    c = c + 0.1;
end

if g_min < 0
    actset = find((gamma > epsilon1) | (g_scale < c * g_min | g < -epsilon2));
else
    actset = find(gamma);
end

end