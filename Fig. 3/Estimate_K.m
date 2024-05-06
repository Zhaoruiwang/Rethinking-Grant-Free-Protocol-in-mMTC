function K = Estimate_K(A,sigma2,sampCov)

s1=A(:,1);
b=s1'*sampCov*s1;
b=real(b);
c=real(s1'*s1);
e=b/c^2;
f=sigma2/c;
K=e-f;
end
