function [y, actset_es,cov_time,iter_time,Coord_total,NumberofGradients] = Random_CD_new(gamma, A, sampCov, sigma2, J,thd,Num_usr)
%{
Coordinate descent algorithm
The termination condition is || [x - \nabla f(x)]_+ - x || < epsilon
%}
[L, N] = size(A);

invSig = eye(L) / sigma2;

k = 0;
epsilon = 1e-3;
N_iter=200;%maximum iteration times

% d_total=zeros(Num_usr,1);%step length of all coordinates
% idx_set=[1:1:Num_usr];%index set for the first iteration

Coord_total=0;%the numbers of total coordinates
NumberofGradients=0;

tic
while k < N_iter
    k = k + 1;

  
    idx_set = randperm(N);

    for iter = 1:N
        idx = idx_set(iter);
        % idx = unidrnd(N);
        a_vec = A(:,idx);

        %==============ML=====
        ainvS = a_vec'*invSig;
        b = ainvS*a_vec;
        b = real(b);
        e = real(ainvS*sampCov*ainvS');
        c = (e - b) / b^2;
        d = max(c, -gamma(idx));
        d=min(d,1-gamma(idx));
%         d_total(idx,k)=d;
        %=============Update===============
        if d ~= 0
            gamma(idx) = gamma(idx) + d;
            g = real(d*b);
            invSig = invSig - d/(1+g) * (ainvS'*ainvS);
        end
    end
%     d_total=[d_total,zeros(Num_usr,1)];
    Coord_total=Coord_total+N;
    gg=gradd(A, sampCov, invSig);
    NumberofGradients=NumberofGradients+N;
    V_gamma=V(N,gamma,gg);
    if norm(V_gamma,'inf') < epsilon
        break
    end
end

iter_time=k;   %times
cov_time=toc;  %time duration

y = gamma;
actset_es=find(gamma);

end




