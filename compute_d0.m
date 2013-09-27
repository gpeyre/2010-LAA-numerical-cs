function d0 = compute_d0(D, x)

% compute_d0 - compute d0 vector
%
%   d0 = compute_d0(D, x);
%
%   d0 =  D1 * ( (D1'*D1)^(-1) * x1 )
%
%   where D1 and x1 are extracted from D and x by selecting non zero
%   entries of x.
%
%   This vector first appears in the work of J.J. Fuchs
%
%       J.J. Fuch, Recovery of exact sparse representations in the presence of bounded noise.
%       IEEE-T-IT, vol.~51, 10, p.~3601--3608, oct.2005. 
%
%   Copyright (c) 2008 Gabriel Peyre

I = find( x~=0 );
D1 = D(:,I);
d0 = D1*( (D1'*D1) \ x(I));