% test for the weak greedy deltaS estimator

path(path, 'toolbox/');

rep = 'results/weak-greedy-deltas-fouclai/';
if not(exist(rep))
    mkdir(rep);
end

n = 100;
p = 2*n;
smax = 20;
methods = {'exact', 'eigen', 'd0'};


n = 1000;
p = 2*n;
smax = 20;
methods = {'eigen', 'd0'};


p = 500;
n = 4*p;
smax = 20;
methods = {'exact', 'd0'};

% A = compute_spherical_matrix(n,p);
A = randn(p,n) / sqrt(p);

%% perform computation
deltas_min = [];
deltas_max = [];
for i=1:length(methods);
    % min
    options.method = ['min' methods{i}];
    disp(['Performing ' options.method '.']);
    deltas_min(:,end+1) = perform_greedy_deltas(A,smax,options);
    % max
    options.method = ['max' methods{i}];
    disp(['Performing ' options.method '.']);
    deltas_max(:,end+1) = perform_greedy_deltas(A,smax,options);
end

%% display deltas
clf;
plot(1:smax, deltas_min);
axis tight;
legend(methods);
saveas(gcf, [rep 'greedy-min-deltas-n' num2str(n) '-p' num2str(p) '.png'], 'png');

clf;
plot(1:smax, deltas_max);
axis tight;
legend(methods);
saveas(gcf, [rep 'greedy-max-deltas-n' num2str(n) '-p' num2str(p) '.png'], 'png');

gr = {'k.-', 'k-', 'k-', 'k*-'};
ms = [15 10 10 10]; 

gr = {'k.', 'k'};
method_lgd = {'Greedy Singular', 'Greedy with d(\sigma)'};
ms = [20 10]; 
lw = 2;

clf;
hold on;
for i=1:2
    h = plot(2:smax,deltas_max(2:end,i), [gr{i} '-']);
    set(h, 'LineWidth', lw);
    set(h, 'MarkerSize', ms(i));
end
for i=1:2
    h = plot(2:smax,deltas_min(2:end,i), [gr{i} '--']);
    set(h, 'LineWidth', lw);
    set(h, 'MarkerSize', ms(i));
end
axis tight; box on;
legend(method_lgd);
saveas(gcf, [rep 'greedy-minmax-deltas-n' num2str(n) '-p' num2str(p) '.eps'], 'eps');

filename = [rep 'greedy-minmax-deltas-n' num2str(n) '-p' num2str(p)];
save(filename, 'deltas_max', 'deltas_min', 'n', 'p', 'method_lgd');