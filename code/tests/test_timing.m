% test for timing

plist = [40 100 200 400 1000];
eta = 4;

plist = [200 1000 5000];
smax = [10 20 40];


plist = [1000 200 500 2000];
smax_list = [30 10 20 40];

options.record_all = 0;
options.verb = 0;
deltas_max = [];
for k=1:length(plist)

    smax = smax_list(k);
    
    p = plist(k);
    n = eta*p;
    A = randn(p,n) / sqrt(p);
    
    disp(['***** n=' num2str(n) ', p=' num2str(p) ', s=' num2str(smax) ' *****']);

    
    if p<=500
    options.extension_size = 1;
    options.pruning_size = n;
    options.method = 'maxexact';
    tic;
    delta = perform_greedy_deltas(A,smax,options);        
    disp(['-> ' options.method ', weak=' num2str(options.extension_size) ...
            ', pruning=' num2str(round(n/options.pruning_size)) ...
            ', deltas=' num2str(delta) ', ' num2str(toc,4) 's.' ]);        
    end
    
    options.extension_size = 1;
    options.pruning_size = n;
    options.method = 'maxd0';
    tic;
    delta = perform_greedy_deltas(A,smax,options);   
    disp(['-> ' options.method ', weak=' num2str(options.extension_size) ...
            ', pruning=' num2str(round(n/options.pruning_size)) ...
            ', deltas=' num2str(delta) ', ' num2str(toc,4) 's.' ]);
        
    options.extension_size = 4;
    options.pruning_size = round(n/10);
    options.method = 'maxd0';
    tic;
    delta = perform_greedy_deltas(A,smax,options);   
    disp(['-> ' options.method ', weak=' num2str(options.extension_size) ...
            ', pruning=' num2str(round(n/options.pruning_size)) ...
            ', deltas=' num2str(delta) ', ' num2str(toc,4) 's.' ]);
        
    options.extension_size = 4;
    options.pruning_size = round(n/100);
    options.method = 'maxd0';
    tic;
    delta = perform_greedy_deltas(A,smax,options);   
    disp(['-> ' options.method ', weak=' num2str(options.extension_size) ...
            ', pruning=' num2str(round(n/options.pruning_size)) ...
            ', deltas=' num2str(delta) ', ' num2str(toc,4) 's.' ]);
       
end