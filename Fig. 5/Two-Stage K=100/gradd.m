function g = gradd(A, sampCov,invSigma)
% Calculate the gradient
% [L, ~] = size(A);
% supp = find(gamma);
% gamma_supp = gamma(supp);
% A_supp = A(:, supp);
% Sigma = (A_supp .* gamma_supp') * A_supp' + sigma2 * eye(L);
% Sigma = (A .* gamma') * A' + sigma2 * eye(L);
% LL = chol(0.5*(Sigma+Sigma'));
% invSigmaA = LL \ (LL' \ A);
invSigmaA =invSigma* A;
B = sampCov * invSigmaA;
C = conj(invSigmaA) .* (A - B);
g = real(sum(C))';
end