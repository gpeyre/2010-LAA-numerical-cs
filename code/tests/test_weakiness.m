% test for the influence of weak-greedy


path(path, 'toolbox/');

rep = 'results/weak-greedy/';
if not(exist(rep))
    mkdir(rep);
end

n = 1000;
p = 2*n;
smax = 20;

A = randn(p,n) / sqrt(p);

weak_list = [1 2 5];


options.pruning_size = n;

deltas_max = [];
deltas_min = [];
lgd = {};
for i=1:length(weak_list)
    q = weak_list(i);
    options.extension_size = q;
    options.method = 'maxd0';
    deltas_max(:,end+1) = perform_greedy_deltas(A,smax,options);    
    options.method = 'mind0';
    deltas_min(:,end+1) = perform_greedy_deltas(A,smax,options);
    lgd{end+1} = ['weak=' num2str(q)];
end


options.extension_size = 1;
options.method = 'maxeigen';
deltas_max(:,end+1) = perform_greedy_deltas(A,smax,options);
options.method = 'mineigen';
deltas_min(:,end+1) = perform_greedy_deltas(A,smax,options);
lgd{end+1} = 'eigen';


clf;
hold on;
plot(deltas_max, '--');
plot(deltas_min, '-');
axis tight;
legend(lgd);
saveas(gcf, [rep 'weak-greedy-minmax-deltas-n' num2str(n) '-p' num2str(p) '.png'], 'png');


% difference

a = length(lgd);
deltas_max_1 = deltas_max(:,1:end-1) - repmat( deltas_max(:,end), [1 a-1] );
deltas_min_1 = deltas_min(:,1:end-1) - repmat( deltas_min(:,end), [1 a-1] );

clf;
plot(deltas_max_1); axis tight;
legend({lgd{1:end-1}});
saveas(gcf, [rep 'weak-greedy-max-deltas-n' num2str(n) '-p' num2str(p) '.png'], 'png');

clf;
plot(deltas_min_1); axis tight;
legend({lgd{1:end-1}});
saveas(gcf, [rep 'weak-greedy-min-deltas-n' num2str(n) '-p' num2str(p) '.png'], 'png');