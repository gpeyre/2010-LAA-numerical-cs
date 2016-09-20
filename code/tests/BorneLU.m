function [VPMAX,VPMIN,mu]=BorneLU(delta,rho)
close all;
pas=0.0001;
t=[pas:pas:1];
t2=[1:pas:3];
temp2=Entr(rho*delta);
%for k=1:length(t)
%  temp(k)=delta*Psimin(t(k),rho)+temp2;
%end
  temp = delta*Psimin(t,rho)+temp2;
%stem(t,temp);
[l,i]=min(abs(temp));
VPMIN=t(i);
%for k=1:length(t2)
%  temp(k)=delta*Psimax(t2(k),rho)+temp2;
%end
  temp =delta*Psimax(t2,rho)+temp2;
%figure;stem(t2,temp,'r');
[l,i]=min(abs(temp));
VPMAX=t2(i);
mu=((1+sqrt(2))/4)*((VPMAX)/(VPMIN)-1);
end
function temp=Psimin(lambda,rho)
temp=Entr(rho)+0.5*((1-rho)*log(lambda)+1-rho+rho*log(rho)-lambda);
end
function temp=Psimax(lambda,rho)
temp=0.5*((1+rho)*log(lambda)+1+rho-rho*log(rho)-lambda);
end
function temp=Entr(p);
temp=-p*log(p)+(p-1)*log(1-p);
end