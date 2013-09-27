function [deltas,deltas_d0] = perform_greedy_deltas(A,smax,options)

% perform_greedy_deltas - greedy deltaS algorithm
%
%   Implement the greedy algorigthm described in the paper 
%
%       A Numerical Exploration of Compressed Sampling Recovery
%       Charles Dossal, Gabriel Peyré and Jalal Fadili,
%       Linear Algebra and its Applications, to appear, 2009.
%       http://hal.archives-ouvertes.fr/hal-00402455/
%
%   deltas = perform_greedy_deltas(A,smax,options);
%
%   A is the (p,n) sensing matrix.
%   deltas(s) is an estimation of the max or min delta_s constant for
%       s<=smax being the sparsity.
%
%   options.method controls wether you want to compute the max or the min
%       delta_s, and which method you would like to use.
%
%   options.method = 'minXXXXX' means you compute the min, while 
%   options.method = 'maxXXXXX' means you compute the max.
%   XXXXX should be replaced by:
%       * XXXXX=exact   -> use a very slow eigenvalue decomposition.
%       * XXXXX=d0      -> use the fast d0 approximation decribed in the paper
%       * XXXXX=eigen   -> use an alternative, slower, method.
%
%   options.pruning_size controls how number of vector are kept at each
%       size. Set it to be much smaller than n if you want speed up (but
%       less accuracy).
%
%   options.extension_size controls how many extension are considered for
%       each vector. Leave this to 1 to achieve fast result, and increase
%       it (e.g. to 3) to have more precision.
%
%   Copyright (c) 2009 Gabriel Peyre, Charles Dossal and Jalal Fadili

options.null = 0;
method = getoptions(options, 'method', 'maxexact');
verb = getoptions(options, 'verb', 1);
record_all = getoptions(options, 'record_all', 1);

n = size(A,2);
p = size(A,1);

% size of the extension
q = getoptions(options, 'extension_size', 1);
% number of kept extension
r = getoptions(options, 'pruning_size', n);
% use sign informartion ?
nosign = getoptions(options, 'nosign', 0);

clb = ['extension_' method];
type = method(1:3);

% initialize supports and signs
Ilist = (1:n)';
Slist = ones(n,1);
deltas = [];
if record_all
    deltas(1) = 0;
end
deltas_d0(1) = 0;
for s=2:smax
    if verb
        progressbar(s-1,smax-1);
    end
    Inew = [];
    Snew = [];
    score = [];
    % greedy extension
    for i=1:size(Slist,1)
        [Iext,Sext,sc] = feval(clb, A, Ilist(i,:), Slist(i,:), q );
        if nosign
            Sext = Sext*0 + 1;
        end
        Inew = [Inew; Iext];
        Snew = [Snew; Sext];
        score = [score; sc(:)];
    end
    % remove doublons
    [Ilist,Slist,score] = unique_support(Inew,Snew,score);
    % Keep only highest score
    [score,sel] = sort(score, 'descend');
    score = score(1:min(end,r));
    sel = sel(1:min(end,r));
    Ilist = Ilist(sel,:); Slist = Slist(sel,:);
    %
    if s==smax || record_all
        % record deltas values
        delta = [];
        for i=1:size(Ilist,1)
            [vmin,vmax] = compute_minmax_eigen(A(:,Ilist(i,:)));
            if strcmp(type, 'max')
                delta(i) = vmax-1;
            else
                delta(i) = 1-vmin;
            end
        end
        [deltas(end+1),i] = max(delta);
    end
    if nargout>1
        % record the bound obtained by d0
        if strcmp(type, 'max')
            deltas_d0(s) = s/norm( my_compute_d0(A(:,Ilist(i,:)),Slist(i,:)) )^2-1;
        else
            deltas_d0(s) = 1 - s/norm( my_compute_d0(A(:,Ilist(i,:)),Slist(i,:)) )^2;
        end
    end
end
deltas = deltas(:);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Iext,Sext,sc] = extension_maxexact(A, I, S, q)
% the larger the score, the better
n = size(A,2);
sc = [];
for i=1:n
    [vmin,sc(i)] = compute_minmax_eigen(A(:,union(I,i)));
end
[sc,i] = max(sc); sc = sc-1;
Iext = [I i];
% compute best sign using previous d0 ...
d0 = my_compute_d0(A(:,I), S);
Sext = [S sign(d0'*A(:,i))];

function [Iext,Sext,sc] = extension_minexact(A, I, S, q)
% the larger the score, the better
n = size(A,2);
sc = [];
for i=1:n
    [sc(i),vmax] = compute_minmax_eigen(A(:,union(I,i)));
end
[sc,i] = min(sc); sc = 1-sc;
Iext = [I i];
% compute best sign using previous d0 ...
d0 = my_compute_d0(A(:,I), S);
Sext = [S -sign(d0'*A(:,i))];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Iext,Sext,sc] = extension_maxd0(A, I, S, q)
s = length(I);
d0 = my_compute_d0(A(:,I),S);
% correlation
c = A'*d0; c(I) = 0;
% best match with sign
% [tmp,i] = min( abs(abs(c)-1) );
[tmp,i] = sort( abs(abs(c)-1), 'ascend');
%
Iext = zeros(q,s+1); Sext = Iext; sc = [];
for k=1:q
    % extend
    Iext(k,:) = [I i(k)];
    Sext(k,:) = [S sign(c(i(k)))];
    % the score is s/|d0|^2-1
    sc(k) = (s+1) / norm( my_compute_d0(A(:,Iext(k,:)),Sext(k,:)) )^2 - 1;
end

function [Iext,Sext,sc] = extension_mind0(A, I, S, q)
s = length(I);
d0 = my_compute_d0(A(:,I),S);
% correlation
c = A'*d0; c(I) = 0;
% best match with sign
% [sc,i] = max(abs(c));
[tmp,i] = sort( abs(c), 'descend');
% Iext = [I i];
% Sext = [S -sign(c(i))];
%
Iext = zeros(q,s+1); Sext = Iext; sc = [];
for k=1:q
    % extend
    Iext(k,:) = [I i(k)];
    Sext(k,:) = [S -sign(c(i(k)))];
    % the score is s/|d0|^2-1
    sc(k) = 1 - (s+1) / norm( my_compute_d0(A(:,Iext(k,:)),Sext(k,:)) )^2;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Iext,Sext,sc] = extension_maxeigen(A, I, S, q)
% compute eigen-matrix
AI = A(:,I);
[U,T] = eig(AI'*AI);
dmax = AI * U(:,end);
% corelation
c = A'*dmax; c(I) = 0;
% extend
[sc,i] = max( abs(c) ); sc = sc-1;
Iext = [I i];
% compute best sign using previous d0 ...
d0 = my_compute_d0(A(:,I), S);
Sext = [S sign(d0'*A(:,i))];

function [Iext,Sext,sc] = extension_mineigen(A, I, S, q)
% compute eigen-matrix
AI = A(:,I);
[U,T] = eig(AI'*AI);
dmin = AI * U(:,1);
% corelation
c = A'*dmin; c(I) = 0;
% extend
[sc,i] = max( abs(c) ); sc = 1-sc;
Iext = [I i];
% compute best sign using previous d0 ...
d0 = my_compute_d0(A(:,I), S);
Sext = [S -sign(d0'*A(:,i))];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d0 = my_compute_d0(AI,S)
Gram=AI'*AI;
d0 = AI*(Gram\S(:));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ilist,Slist,score] = unique_support(Ilist,Slist,score)
% sort all support
for i=1:size(Ilist,1)
    [Ilist(i,:),J] = sort(Ilist(i,:),2);
    Slist(i,:) = Slist(i,J);
end
% select unique support
[Ilist,J] = unique(Ilist,'rows');
Slist = Slist(J,:);
score = score(J);

%%

function [Iext,Sext,sc] = extension_maxd0_charles(A, I, S, q)
nbcandidats = 1;
d1=d02(A,I,S);
ps=A'*d1;
ps(I)=0;
[pss,Ordre]=sort(abs(ps),'descend');
Ind2=Ordre(1:nbcandidats);
si2=sign(ps(Ind2));

Iext = [I, Ind2(1)];
Sext = [S, si2(1)];
sc = 0;


function Res=d02(A,Ind,si)
AI=A(:,Ind);
Gram=AI'*AI;
Res=AI*(Gram\si');