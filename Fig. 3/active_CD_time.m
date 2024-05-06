function [y, actset_es,cov_time,iter_time] = active_CD_time(gamma, A, sampCov, sigma2, J,Num_usr, max_time)

[L, N] = size(A);

invSig = eye(L) / sigma2;
k = 0;
epsilon = 1e-3;
allset=1:1:N;
cov = tic;
while toc(cov) <=max_time
    k = k + 1;

    if k==1
        V_gamma=V_pick(A,sampCov,invSig,gamma,allset);
        omega=max(0.2^(k+1)*norm(V_gamma,'inf'),epsilon);
        active_set=find(V_gamma>=omega);
    else
        omega=max(0.2^(k+1)*norm(V_gamma,'inf'),epsilon);
        active_set=find(V_gamma>=omega);
    end
    %% iteration

    for  indx = active_set
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
    V_gamma=V_pick(A,sampCov,invSig,gamma,allset);

    % if norm(V_gamma,'inf') < epsilon
    %     break
    % end
end

iter_time=k;   %times
cov_time=toc(cov);  %time duration

y = gamma;
actset_es=find(gamma);

end
