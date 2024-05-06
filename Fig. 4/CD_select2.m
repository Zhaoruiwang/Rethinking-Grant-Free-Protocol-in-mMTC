function [y, actset_es,cov_time,iter_time,Coord_total,NumberofGradients] = CD_select2(gamma, A, sampCov, sigma2, J,thd,Num_usr,K_est,threshold,count_time)
%{
Coordinate descent algorithm
The termination condition is || [x - \nabla f(x)]_+ - x || < epsilon
%}
[L, N] = size(A);
invSig = eye(L) / sigma2;
k = 0;
epsilon = 1e-3;
N_iter=200;%maximum iteration times
Coord_total=0;
NumberofGradients=0;
cov=tic;
count=zeros(N,1);
while k < N_iter
    %     startloop=tic;
    k = k + 1;
    %% select
    if k==1
        g_act=gradd(A,sampCov,invSig);
        NumberofGradients=NumberofGradients+N;
        V_gamma=V(N,gamma,g_act);
        [~, sorted_indices]=sort(V_gamma,'descend');
        selected_set=sorted_indices(1:K_est);
    else
        [~, sorted_indices]=sort(V_gamma,'descend');
        if length(idx_set)>=K_est
            selected_set=sorted_indices(1:K_est);
        else
            selected_set=sorted_indices(1:length(idx_set));
        end
    end

    %% iteration
    for  indx = selected_set
        a_vec = A(:,indx);
        %==============ML=====
        ainvS = a_vec'*invSig;
        b = ainvS*a_vec;
        b = real(b);
        e = real(ainvS*sampCov*ainvS');
        c = (e - b) / b^2;
        % c = real(ainvS*sampCov*ainvS' - b) / abs(b).^2;
        d = max(c, -gamma(indx));
        d=min(d,1-gamma(indx));
        %=============Update===============
        if d ~= 0
            gamma(indx) = gamma(indx) + d;
            g = real(d*b);
            invSig = invSig - d/(1+g) * (ainvS'*ainvS);
        end
    end
    Coord_total=Coord_total+length(selected_set);
    zero_idx=find(V_gamma<=threshold);
    count(zero_idx)=count(zero_idx)+1;
    idx_set=find(count<=count_time);
    V_gamma=V_pick(A,sampCov,invSig,gamma,idx_set);
    NumberofGradients=NumberofGradients+length(idx_set);
    if norm(V_gamma,'inf') < epsilon
        break
    end
end
iter_time=k;   %times
cov_time=toc(cov);  %time dur
y=gamma;
actset_es=find(gamma);
end