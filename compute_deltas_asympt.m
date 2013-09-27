function [delta_max, delta_min, mu]=compute_deltas_asympt(a,b,c)

%   compute_deltas_asympt - compute the asymptotics of delta_s.
%
%   [delta_max, delta_min, mu]=compute_deltas_asympt(delta,rho);
% or
%   [delta_max, delta_min, mu]=compute_deltas_asympt(s,p,n)
% with s<p<n
%   delta=p/n, rho=s/p
%
%   Compute the lower bounds found on the paper
%       Jeffrey D. Blanchard, Jarred Tanner and Coralia Cartis 
%       Decay properties for restricted isometry constants;  
%       IEEE Signal Processing Letters, Vol. 16(7) (2009) 572-575
%
%   Copyright (c) 2009 Charles Dossal and Gabriel Peyre

if nargin==2
    delta = a;
    rho = b;
else
    s = a;
    p = b;
    n = c;
    rho = s/p;
    delta = p/n;
end


if length(rho)>1
    for i=1:length(rho)
        [delta_max(i), delta_min(i), mu(i)]=compute_deltas_asympt(delta,rho(i));
    end
    return;
end


nstep = 10000;
t = linspace(0,1,nstep); t(1)=[];
t2 = linspace(1,3,nstep);

temp2=Entr(rho*delta);
temp = delta*Psimin(t,rho)+temp2;
[l,i]=min(abs(temp));
VPMIN=t(i);
temp = delta*Psimax(t2,rho)+temp2;
[l,i]=min(abs(temp));
VPMAX=t2(i);
mu=((1+sqrt(2))/4)*((VPMAX)/(VPMIN)-1);

delta_max = VPMAX - 1;
delta_min = 1 - VPMIN;

end

function temp=Psimin(lambda,rho)
temp=Entr(rho)+0.5*((1-rho)*log(lambda)+1-rho+rho*log(rho)-lambda);
end

function temp=Psimax(lambda,rho)
temp=0.5*((1+rho)*log(lambda)+1+rho-rho*log(rho)-lambda);
end

function temp=Entr(p);
temp=-p.*log(p)+(p-1).*log(1-p);
end