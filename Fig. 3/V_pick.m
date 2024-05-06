function V_gamma = V_pick(A,sampCov,invSigma,gamma,idx_set)
    N=length(gamma);


    gg=gradd(A(:,idx_set), sampCov,invSigma);
    V_gamma=V1(N,idx_set,gamma(idx_set),gg);
end