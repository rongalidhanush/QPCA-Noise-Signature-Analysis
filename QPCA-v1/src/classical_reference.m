function out = classical_reference(rho)
%CLASSICAL_REFERENCE Numerical reference for evaluating QPCA/QPE.
%
% This function is NOT the quantum algorithm. It provides ground-truth
% eigenvalues/eigenvectors for validation.

rho = (rho+rho')/2;
[V,D] = eig(rho);
lambda = real(diag(D));
[lambda,idx] = sort(lambda,'descend');
V = V(:,idx);

lambda(lambda<0)=0;
lambda=lambda/sum(lambda);

out.eigenvalues=lambda;
out.eigenvectors=V;
out.explainedPercent=100*lambda;
out.cumulativePercent=cumsum(100*lambda);

end
