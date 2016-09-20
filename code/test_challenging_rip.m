% test the RIP estimation by the greedy algorithm.

%%
% Add the helpers to the path.

path(path, 'toolbox/');

%% 
% Parameters for the problems.

% number of measures
p = 500;
% dimension of the problem
n = 4*p;
% maximum tested sparsity
smax = 12;
% random sensing matrix.
A = randn(p,n) / sqrt(p);

%%
% To speed up the algorithm, we use a 1/8 pruning. To increase its
% precision, we use a x2 extension

options.extension_size = 2;
options.pruning_size = round(n/8);

%%
% Compute minimum RIP constant using the greedy algorithm.

options.method = 'mind0';
deltas_min = perform_greedy_deltas(A,smax,options);


%%
% Compute maximum RIP constant using the greedy algorithm.

options.method = 'maxd0';
deltas_max = perform_greedy_deltas(A,smax,options);

%%
% Compute theoritical upper bounds

[deltas_max_th, deltas_min_th, mu]=compute_deltas_asympt(1:smax,p,n);

%% 
% Display curves.
clf;
hold on;
plot(1:smax, deltas_min, 'k-');
plot(1:smax, deltas_max, 'k.-');
plot(1:smax, deltas_min_th, 'k--');
plot(1:smax, deltas_max_th, 'k.--');
axis tight;
legend('Min', 'Max', 'Min asymp', 'Max asymp');
