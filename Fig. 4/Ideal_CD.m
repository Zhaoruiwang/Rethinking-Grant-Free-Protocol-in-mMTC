function [y, act_set_ideal,cov_time, iter_time,Coord_total] = Ideal_CD(gamma,A, sampCov, sigma2, actset,Num_usr)
[L, N] = size(A);

invSig = eye(L) / sigma2;

A_act = A(:, actset);
gamma_act = gamma(actset);
N=[1:1:Num_usr];
Coord_total=0;%the numbers of total coordinates

r = length(actset);


epsilon = 1e-3;
N_iter=200;

tic;
k = 0;

while k < N_iter
    k = k + 1;
    
    idx_set = randperm(r);
    for iter = 1:r
        idx = idx_set(iter);
        % idx = unidrnd(N);
        a_vec = A_act(:,idx);
        
        %==============ML=====
        ainvS = a_vec'*invSig;
        b = ainvS*a_vec;
        b = real(b);
        e = real(ainvS*sampCov*ainvS');
        c = (e - b) / b^2;
        % c = real(ainvS*sampCov*ainvS' - b) / abs(b).^2;
        d = max(c, -gamma_act(idx));
        %=============Update===============
        if d ~= 0
            gamma_act(idx) = gamma_act(idx) + d;
            g = real(d*b);
            invSig = invSig - d/(1+g) * (ainvS'*ainvS);
        end
    end
    gamma(actset) = gamma_act;
    Coord_total=Coord_total+r;
    g=grad(gamma, A, sampCov, sigma2);
    if norm(max(gamma - g, 0) - gamma) < epsilon
        break
    end
end

cov_time=toc;
iter_time=k;

gamma(actset) = gamma_act;
y = gamma;
act_set_ideal=find(gamma);
end