function [vmin,vmax] = compute_minmax_eigen(A)

% compute_minmax_eigen - compute the min/max eigenvalues of A'*A
%
%   [vmin,vmax] = compute_minmax_eigen(A);
%
%   Copyright (c) 2008 Gabriel Peyre

[U,S] = eig(A'*A);
vmin = min(diag(S));
vmax = max(diag(S));

return;

[U,S,V] = svd(A);
vminX(end+1) = min(diag(S))^2;
vmaxX(end+1) = max(diag(S))^2;