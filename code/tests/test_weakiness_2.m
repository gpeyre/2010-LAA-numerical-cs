% test for the influence of weak-greedy

path(path, 'toolbox/');

rep = 'results/weak-greedy/';
if not(exist(rep))
    mkdir(rep);
end

p = 1000;
n = 4*p;
smax = 20;

p = 3000;
n = 4*p;
smax = 30;


p = 2000;
n = 4*p;
smax = 30;

A = randn(p,n) / sqrt(p);

weak_list = [1 4];
prune_list = [n round(n/4)];


deltas_max = [];
deltas_min = [];
method_lgd = {};
for i=1:length(weak_list)
    r = weak_list(i);
    q = prune_list(i);
    options.extension_size = r;
    options.pruning_size = q;
    options.method = 'maxd0';
    deltas_max(:,end+1) = perform_greedy_deltas(A,smax,options);    
    options.method = 'mind0';
    deltas_min(:,end+1) = perform_greedy_deltas(A,smax,options);
    method_lgd{end+1} = ['N/Q=' num2str(n/q), ', R=' num2str(r)];
end


smaxd = smax;

gr = {'k.', 'k*', 'k'};
ms = [20 10 15]; 
lw = 2;

[deltas_max(:,end+1), deltas_min(:,end+1)]=compute_deltas_asympt(1:smax,p,n);
method_lgd{end+1} = 'Asymptotic upper bound';

clf;
hold on;
for i=1:3
    h = plot(2:smaxd,deltas_max(2:smaxd,i), [gr{i} '-']);
    set(h, 'LineWidth', lw);
    set(h, 'MarkerSize', ms(i));
end
for i=1:3
    h = plot(2:smaxd,deltas_min(2:smaxd,i), [gr{i} '--']);
    set(h, 'LineWidth', lw);
    set(h, 'MarkerSize', ms(i));
end
    
axis tight; box on;
legend(method_lgd);
saveas(gcf, [rep 'wgreedy-minmax-n' num2str(n) '-p' num2str(p) '.eps'], 'eps');
saveas(gcf, [rep 'wgreedy-minmax-n' num2str(n) '-p' num2str(p) '.png'], 'png');

%% display Foucart/Lai bound, delta_S style
% (4\sqrt{2}-3) \ripl{2s} + \ripu{2s} < 4 (\sqrt{2}-1),
slist = 2:smaxd;
clf;
hold on;
for i=1:3
    L = deltas_min(slist,i);
    U = deltas_max(slist,i);    
    T = (1+sqrt(2))*L + U;
    T = (4*sqrt(2)-3) * L + U;
    h = plot( slist/2, T, [gr{i} '-'] );
    set(h, 'LineWidth', lw);
    set(h, 'MarkerSize', ms(i));
end
legend(method_lgd);
h = plot([slist(1) slist(end)]/2, 4*(sqrt(2)-1)*[1 1], 'k--');
set(h, 'LineWidth', 2);
hold off;
axis tight; box on;
saveas(gcf, [rep 'wgreedy-fouclai-n' num2str(n) '-p' num2str(p) '.eps'], 'eps');
saveas(gcf, [rep 'wgreedy-fouclai-n' num2str(n) '-p' num2str(p) '.png'], 'png');

