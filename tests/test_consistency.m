% test for the consistency of the estimation

path(path, 'toolbox/');

rep = 'results/consistency/';
if not(exist(rep))
    mkdir(rep);
end

% delta=p/n, rho=s/p
rho = 1/100;

ntrials = 20;

slist = 3:8;
delta_list = [1/2 1/4];

deltas_min = zeros(length(slist), length(delta_list));
deltas_max = deltas_min;

for i=1:length(slist)
    %    disp(['delta=' num2str(delta)]);
    for j=1:length(delta_list)
        for k=1:ntrials
            delta = delta_list(j);
            s = slist(i);
            p = round( s/rho );
            n = round( p/delta );
            A = randn(p,n) / sqrt(p);
            %
            options.extension_size = 4;
            options.pruning_size = round(n/4);
            options.record_all = 0;
            %
            options.method = 'maxd0';
            deltas_max(i,j) = deltas_max(i,j) + perform_greedy_deltas(A,s,options);
            options.method = 'mind0';
            deltas_min(i,j) = deltas_min(i,j) + perform_greedy_deltas(A,s,options);
        end
    end
end
deltas_max = deltas_max/ntrials;
deltas_min = deltas_min/ntrials;

lgd = {'\eta=1/2' '\eta=1/4'};
gr = {'k.', 'k*', 'k'};
ms = [20 10 15]; 
lw = 2;

clf;
hold on;
for i=1:length(delta_list)
    h = plot(slist,deltas_max(:,i), [gr{i} '-']);
    set(h, 'LineWidth', lw);
    set(h, 'MarkerSize', ms(i));
end
for i=1:length(delta_list)
    h = plot(slist,deltas_min(:,i), [gr{i} '--']);
    set(h, 'LineWidth', lw);
    set(h, 'MarkerSize', ms(i));
end    
axis tight; box on;
legend(lgd);
saveas(gcf, [rep 'consistency-n' num2str(n) '-p' num2str(p) '.eps'], 'eps');
saveas(gcf, [rep 'consistency-n' num2str(n) '-p' num2str(p) '.png'], 'png');
