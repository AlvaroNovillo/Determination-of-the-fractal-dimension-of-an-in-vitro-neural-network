%WattsStrogatz(N,K,beta)
%N numero de nodos
%K average degree: número de enlaces que tiene con otros nodos.
% beta es la probabilidad de "recablear" los enlaces inicialmente asignados a un nodo en la estructura anillo
% beta = 0 is a ring lattice, and beta = 1 is a random graph.
%h2 = WattsStrogatz(10,1,0)
openExample('matlab/BuildWattsStrogatzSmallWorldGraphModelExample')
%% 
% Study of C. elegans network using W-S model.

clear all;
SmallWorld(282,14)
%%
%Watts Strogatz rewiring method visualization
h2 = WattsStrogatz(50,3,0.0000000);
figure();
hold on;
t = tiledlayout(1,3,'TileSpacing','Compact');
nexttile
plot(h2,'o','MarkerSize',4,'NodeLabel', {},'NodeColor','k','EdgeColor',[17 17 17]/255);
title('$p = 0$','Interpreter','latex')
nexttile
h3 = WattsStrogatz(50,3,10^(-1.5));
plot(h3,'o','MarkerSize',4,'NodeLabel', {},'NodeColor','k','EdgeColor',[17 17 17]/255);
title('$p = 0.0316$','Interpreter','latex')
nexttile
h4 = WattsStrogatz(50,3,0.5);
plot(h4,'o','MarkerSize',4,'NodeLabel', {},'NodeColor','k','EdgeColor',[17 17 17]/255);
hold off;
title('$p = 0.5$','Interpreter','latex')
title(t, 'Watts-Strogatz Graph with $N = 50$ nodes, $K = 3$', ...
    'Interpreter','latex')
    
%% 
% 

function SmallWorld(n,k)
%Numerical simulation of Watts-Strogatz model
%Calculamos el factor de normalización

G_0 = WattsStrogatz(n,k,0);
A_0 = adjacency(G_0);
d_0 = distances(G_0);
C_0 = clustering_coef_bu(A_0);
d_mean_0 = mean(mean(d_0));
C_mean_0 = mean(C_0);

%Creamos un vector p que esta espaciado logaritmicamente
%(para enfatizar el cambio del clustering y del shortest path length en valores muy próximos...
% ... a beta=0)

p=logspace(1.0022,10,31).*10^(-10);
colToDelete = p < 0.0001 ;
p(colToDelete) = [];

%Calculamos L/L(0) y C/C(0) para cada valor de beta:

C_mean_array = [];
d_mean_array = [];
C_total = [];
d_total = [];

for j=1:20

    for i = 1:length(p) 
        G = WattsStrogatz(n,k,p(1,i));
        A = adjacency(G);
        d = distances(G);
        C = clustering_coef_bu(A);
        d_mean = mean(mean(d))/d_mean_0;
        C_mean=mean(C)/C_mean_0;
        C_mean_array(j,i) = C_mean;
        d_mean_array(j,i) = d_mean;
        i=i+1;
    end

end
C_mean_total = mean(C_mean_array);
d_mean_total = mean(d_mean_array);


%Ploteamos los resultados, usando escala logarítmica en el eje de las x:

figure();
hold on;
plot(p(1,:),d_mean_total(1,:),'o')
plot(p(1,:), C_mean_total(1,:),'s')
set(gca,'Xscale','log')
yticks([0,0.2,0.4,0.6,0.8,1])
legend('$L(p)/L(0)$','$C(p)/C(0)$','Interpreter','latex')
hold off 

end

function C=clustering_coef_bu(G)

%CLUSTERING_COEF_BU     Clustering coefficient
%
%   C = clustering_coef_bu(A);
%
%   The clustering coefficient is the fraction of triangles around a node
%   (equiv. the fraction of node?s neighbors that neighbors of each other).
%
%   Input:      A,      binary undirected connection matrix
%
%   Output:     C,      clustering coefficient vector
%
%   Reference: Watts and Strogatz (1998) Nature 393:440-442.
%
%
%   Mika Rubinov, UNSW, 2007-2010
n=length(G);
C=zeros(n,1);
for u=1:n
    V=find(G(u,:));
    k=length(V);
    if k>=2                %degree must be at least 2
        S=G(V,V);
        C(u)=sum(S(:))/(k^2-k);
    end
end
end