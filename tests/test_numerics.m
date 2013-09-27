% test for computation speed of eigenvectors

slist = round(linspace(3,100,20));
p = 200;
ntrial = 100;

for i=1:length(slist)
    progressbar(i,length(slist));
    s = slist(i);
    AI = randn(p,s)/sqrt(p);
    S = sign(randn(s,1));
    % computation of d0
    tic;
    for k=1:ntrial
        G = AI'*AI;
        d0 = AI*(G\S);
    end
    Td0(i) = toc/ntrial;
    % eigenvalue max
    tic;
    for k=1:ntrial
        G = AI'*AI;
        a = svds(G, 1, 'L');
    end
    TsvdsMax(i) = toc/ntrial;    
    % eigenvalue min
    tic;
    for k=1:ntrial
        G = AI'*AI;
        a = svds(G, 1, 0);
    end
    TsvdsMin(i) = toc/ntrial;  
    % eigenvalue minmax
    tic;
    for k=1:ntrial
        G = AI'*AI;
        [S,D] = eig(G);
    end
    Teig(i) = toc/ntrial;    
end

clf;
subplot(2,1,1);
hold on
plot(slist, TsvdsMax ./ Td0, 'b');
plot(slist, TsvdsMin ./ Td0, 'r');
axis tight;
legend('SVDS, max', 'SVDS, min');
subplot(2,1,2);
plot(slist, Teig ./ Td0, 'g');
legend('Eig');

rep = 'results/numerics/';
if not(exist(rep))
    mkdir(rep);
end
saveas(gcf, [rep 'numerics.png'], 'png');
