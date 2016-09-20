% test for bug
load MatriceProblematique

% n = size(A,2);
% p = size(A,1);

p = 400;
n = 2*p;

A = randn(p,n)/sqrt(p);
smax = 8;

nbcandidats = 1;
nbVecteurs = n;
[VPMAX,Approxd0]=BattreBP8(A,smax,nbVecteurs,nbcandidats);


options.method = 'maxd0';
[deltas,deltas_d0] = perform_greedy_deltas(A,smax,options);

VPMAX-1
deltas(end)