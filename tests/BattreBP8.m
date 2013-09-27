function [VPMAX,Approxd0]=BattreBP8(A,d,nbVecteurs,nbcandidats)
% Le programme BattreBP renvoie des vecteurs non identifiables associés à
% une matrice A qui ont une norme l_0 égale à d.
% A est une matrice de taille n times p.
% d est un entier représentant la sparsité à atteindre.
% nbVecteurs représente le nombre de Vecteurs étudiés par l'algorithme à chaque étape.
% Une étape correpsondant à une sparsité k<=d.
% En pratique on peut prendre 50 pour commencer mais on peut monter à 400
% sans trop de pb.
% nbcandidats représente le nombre de manière dont est étendu chaque
% vecteur courant. En pratique nbcandidats varie entre 5 et 10.
% VPmax est la VP maximale d'une sous matrice associée.
[n,p]=size(A);
%I=[1:p]';
I=[1:p]';
%SI=ones(p,1);
SI=ones(p,1);
for i=1:d-1
    %On r?ordonne les Supports et les Signes associ?s
    %pour pouvoir supprimer les doublons.
    n1=size(I,1);
    [I2,Ind]=sort(I,2);
    for j=1:n1
        SI2(j,:)=SI(j,Ind(j,:));
    end
    I=I2;
    SI=SI2;
    %Suppression des doublons
    [I,J]=unique(I,'rows');
    SI=SI(J,:);
    %On ne garde que les vecteurs ayant une grande norme de d0.
    n1=size(I,1);
    D0=zeros(n,n1);
    NormeD0=zeros(n1,1);
    for j=1:n1
        D0(:,j)=d02(A,I(j,:),SI(j,:));
        NormeD0(j)=norm(D0(:,j),2);
    end
    if (n1>nbVecteurs)
        [temp,Ind]=sort(NormeD0,'ascend');
        temp=temp(1:nbVecteurs);
        Ind=Ind(1:nbVecteurs);
        D0=D0(:,Ind);
        I=I(Ind,:);
        SI=SI(Ind,:);
    end
    %Calcul des nouveaux supports et signes.
    n1=size(I,1);
    I2=[];
    SI2=[];
    for j=1:n1
        Ind=I(j,:);
        si=SI(j,:);
        d1=D0(:,j);
        ps=A'*d1;
        ps(Ind)=0;
        [pss,Ordre]=sort(abs(ps),'descend');
        Ind2=Ordre(1:nbcandidats);
        si2=sign(ps(Ind2));
        for jj=1:nbcandidats
            I2=[I2;I(j,:),Ind2(jj)];
            SI2=[SI2;SI(j,:),si2(jj)];
        end
    end
    I=I2;
    SI=SI2;
end

%Maintenant on s?lectionne les meilleurs donc
%on reprend la proc?dure de trie.
n1=size(I,1);
[I2,Ind]=sort(I,2);
for j=1:n1
    SI2(j,:)=SI(j,Ind(j,:));
end
I=I2;
SI=SI2;
%Suppression des doublons
[I,J]=unique(I,'rows');
SI=SI(J,:);
%On ne garde que les vecteurs ayant une grande norme de d0.
n1=size(I,1);
D0=zeros(n,n1);
NormeD0=zeros(n1,1);
for j=1:n1
    D0(:,j)=d02(A,I(j,:),SI(j,:));
    NormeD0(j)=norm(D0(:,j),2);
end
if (n1>nbVecteurs)
    [temp,Ind]=sort(NormeD0,'ascend');
    temp=temp(1:nbVecteurs);
    Ind=Ind(1:nbVecteurs);
    D0=D0(:,Ind);
    I=I(Ind,:);
    SI=SI(Ind,:);
end



n1=size(I,1);
for k=1:n1
    temp(k)=norm(d02(A,I(k,:),SI(k,:)),2);
end
% stem(temp);hold on;
temp2=0*temp;
x=[];
%h=waitbar(0,'Please Wait ...');
nbtests=n1;
VPMAX=0;
Approxd0=0;
for k=1:nbtests
    AI=A(:,I(k,:));
    temp3=svds(AI'*AI,1);
    if temp3>VPMAX
        VPMAX=temp3;
        Approxd0=d/(temp(k))^2;
    end
end

% close(h);
end

function Res=d02(A,Ind,si);
AI=A(:,Ind);
Gram=AI'*AI;
Res=AI*(Gram\si');
end