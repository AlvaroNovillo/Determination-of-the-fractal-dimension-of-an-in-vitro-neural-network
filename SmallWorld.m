%WattsStrogatz(N,K,beta)
%N numero de nodos
%K average degree: número de enlaces que tiene con otros nodos.
% beta es la probabilidad de "recablear" los enlaces inicialmente asignados a un nodo en la estructura anillo
% beta = 0 is a ring lattice, and beta = 1 is a random graph.
%%
SmallWorld(282,14)

h2 = WattsStrogatz(282,14,0.15);
plot(h2,'NodeColor','k','EdgeAlpha',0.2);
title('Watts-Strogatz Graph with $N = 282$ nodes, $K = 14$, and $\beta = 0.15$', ...
    'Interpreter','latex')
    
%%
function SmallWorld(n,k)

openExample('matlab/BuildWattsStrogatzSmallWorldGraphModelExample')

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

hold on
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