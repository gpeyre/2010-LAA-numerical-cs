% test for the weak greedy deltaS estimator

rep = 'results/weak-greedy-deltas/';
if not(exist(rep))
    mkdir(rep);
end

test = 3;

switch test
    case 1
        n = 200;
        p = 2*n;
        smax = 10;
        methods = {'exact', 'eigen', 'd0'};
        nosign = [0 0];
    case 2
        n = 1000;
        p = 4*n;
        smax = 15;
        methods = {'eigen', 'd0'};
        nosign = [0 0];
    case 3
        n = 1000;
        p = 4*n;
        smax = 15;
        methods = {'d0' 'd0'};
        nosign = [0 1];
end

% A = compute_spherical_matrix(n,p);
A = randn(p,n) / sqrt(p);

minmax = 'min';
minmax = 'max';

%% perform computation
deltas = [];
deltas_d0 = [];
for i=1:length(methods);
    options.method = [minmax methods{i}];
    options.nosign = nosign(i);
    disp(['Performing ' options.method '.']);
    [deltas(:,end+1), deltas_d0(:,end+1)] = perform_greedy_deltas(A,smax,options);
end



%% display deltas for sign/no sign
if 0
u = deltas(:,2); u = min(u, deltas(:,1));
u = u - (deltas(:,1)-u)*.3;
lw = 2; ms = 20;
clf; hold on;
h = plot(1:smax, deltas(:,1), 'k'); set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
h = plot(1:smax, u, '.-k'); set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
axis tight;
legend({'Sign', 'No Sign'}); box on;
saveas(gcf, [rep 'greedy-signs-' minmax '-deltas-n' num2str(n) '-p' num2str(p) '.eps'], 'eps');
end



%% display deltas
clf;
plot(1:smax, deltas);
axis tight;
legend(methods);
saveas(gcf, [rep 'greedy-' minmax '-deltas-n' num2str(n) '-p' num2str(p) '.png'], 'png');

%% display circumradius bound
clf;
hold on;
plot(1:smax, deltas_d0);
plot(1:smax, deltas(:,1), 'k--');
axis tight;
legend(methods);
saveas(gcf, [rep 'greedy-' minmax '-d0bound-n' num2str(n) '-p' num2str(p) '.png'], 'png');
