function [y, actset_es,cov_time,iter_time,d_cd] = Random_CD_Niter(gamma, A, sampCov, sigma2, J,thd,Num_usr,N_iter)
%{
Coordinate descent algorithm
The termination condition is || [x - \nabla f(x)]_+ - x || < epsilon
%}
[L, N] = size(A);


invSig = eye(L) / sigma2;


k = 0;

epsilon = 1e-3;
fvalue_now = f(gamma, A, sampCov, sigma2);

d_total=zeros(Num_usr,1);
idx_set=[1:1:Num_usr];

tic
while k < N_iter
    k = k + 1;
    fvalue_last=fvalue_now;
    
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
        % c = real(ainvS*sampCov*ainvS' - b) / abs(b).^2;
%         if k >1 && gamma(idx)==0
%             idx
%             pause(10);
%         end
        d = max(c, -gamma(idx));
        d_total(idx,k)=d;
        %=============Update===============
        if d ~= 0
            gamma(idx) = gamma(idx) + d;
            g = real(d*b);
            invSig = invSig - d/(1+g) * (ainvS'*ainvS);
        end
    end
    fvalue_now=f(gamma, A, sampCov, sigma2);
    abs_value=abs(fvalue_last-fvalue_now);
    d_total=[d_total,zeros(Num_usr,1)];
    if abs_value < epsilon
        break
    end
end

iter_time=k;   %times
cov_time=toc;  %time duration
% Blk_lth=Q_max+1;
% Blk_num=N/(Q_max+1);
% for b=1:Blk_num
%     gamma_b=gamma((b-1)*Blk_lth+1:b*Blk_lth);
%     [gamma_b_max index]=max(gamma_b);
%     if gamma_b_max>=thd
%         set=(b-1)*Blk_lth+1:b*Blk_lth;
%         set(index)=[];
%         gamma(set)=0;
%     else
%         gamma((b-1)*Blk_lth+1:b*Blk_lth)=0;
%     end
% end

d_cd=d_total;
y = gamma;
actset_es=find(gamma);

end




