clear all
DIVS = [0	1	2	3	4	5	6	7	8	9	10	11	12	13	14 16 18];
%r= rango scaling
for i=1:length(DIVS)
    dum=sprintf('dataPLOS/*DIV%d_*.mat',DIVS(i));
    files=dir(dum);
    nf=length(files);
    for n=1:nf
        filename=horzcat(files(n).folder,'/',files(n).name);
        data=load(filename);
        
        %extraes la red
        AF=data.net.FULL_ADJACENCY; %Matriz de adyacencia completa, que incluye neuronas + bifurcaciones
        AC=data.net.CLUSTER_ADJACENCY;%Matriz de adyacencia "efectiva" que solo incluye neuronas y conexiones directas e indirectas a traves de un pto de bifurcacion
        %Posición espacial de los nodos 
        xc=data.net.CLUSTER_CENTROID(:,1);
        yc=data.net.CLUSTER_CENTROID(:,2);
        
        %Buscamos el subgrafo con las componentes conexas  
        GF=graph(AF,'omitselfloops');
        %El subgrafo de la GF
        [bin,binsize] = conncomp(GF);
        idx = binsize(bin) == max(binsize);
        GCF = subgraph(GF, idx);
                    
        %Adjacency matrix for the Full Graph
        AdGCF = adjacency(GCF); %De la componente gigante
                  
        GC=graph(AC, 'omitselfloops');
        GF=graph(AF,'omitselfloops');

        %El subgrafo de la GC
        [bin,binsize] = conncomp(GC);
        idx = binsize(bin) == max(binsize);
        GCC = subgraph(GC, idx); 
        
        %El subgrafo de la GF
        [bin,binsize] = conncomp(GF);
        idx = binsize(bin) == max(binsize);
        GCF = subgraph(GF, idx);
        
        %Normalizamos por el numero de nodos:
        N_GC = numnodes(GCC);
        N_GF = numnodes(GCF);

        AdGCC = adjacency(GC); %Extendido a todo el grafo
        D_GC= distances(GCC,'Method','unweighted');
        %Sustituimos los valores inf por la media
        D_GC(D_GC==Inf)=0;
        dist_GC(i)= mean(D_GC(D_GC>0))/(N_GC-1);
        C_GC(i) = mean(clustering_coef_bu(AdGCC));

        %Notamos que la normalización es N-1 (Todos los nodos menos el
        %que estas contando)

        %Para la GF
        AdGCF = adjacency(GF); %Extendido a todo el grafo
        D_GF= distances(GCF,'Method','unweighted');
        %Sustituimos los valores inf por la media
        D_GF(D_GF==Inf)=0;
        dist_GF(i)= mean(D_GF(D_GF>0))/(N_GF-1);
        C_GF(i) = mean(clustering_coef_bu(AdGCF));

    end
    errC_GC(i) = std(C_GC)./sqrt(length(C_GC));
    C_GC_mean(i) = mean(C_GC);
    
    
    errC_GF(i) = std(C_GF)./sqrt(length(C_GF));
    C_GF_mean(i) = mean(C_GF);
    
    errdist_GC(i) = std(dist_GC)./sqrt(length(dist_GC));
    dist_GC_mean(i) = mean(dist_GC);
    
    
    errdist_GF(i) = std(dist_GF)./sqrt(length(dist_GF));
    dist_GF_mean(i) = mean(dist_GF);
end
%%
topo_paint(DIVS,C_GC_mean,dist_GC_mean,errC_GC,errdist_GC,C_GF_mean,dist_GF_mean,errC_GF,errdist_GF)
%%
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

function topo_paint(DIVS,C_GC_mean,dist_GC_mean,errC_GC,errdist_GC,C_GF_mean,dist_GF_mean,errC_GF,errdist_GF)
    figure; hold on;
    % GC Mean cluster and shortest path
    [ax,h1,h2]=plotyy(DIVS,C_GC_mean,DIVS,dist_GC_mean,'plot');
    
    set(h1,'marker','o');
    set(h2,'marker','s');
    
    legend('C','L/N','location','best');
    
    hold(ax(1),'on');
    
    e1 = errorbar(ax(1),DIVS,C_GC_mean,errC_GC,'o');
    e1.Color = [0 0.4470 0.7410];
    set(ax(1),'YColor',[0 0 0])
%     set(ax(1),'xlim',[0 max(DIVS)],'ylim', [0,0.3],'yticklabel',[0:0.3],'ytick',[0:0.3])
    ylabel(ax(1),'Mean Clustering Coeff.')
    
    hold off;
    
    hold(ax(2),'on');
    
    e2 = errorbar(ax(2),DIVS,dist_GC_mean,errdist_GC,'s');
    e2.Color = [0.8500 0.3250 0.0980];
    set(ax(2),'YColor',[0 0 0])
%     set(ax(2),'xlim',[0 max(DIVS)],'ylim', [0,0.3],'yticklabel',[0:0.3],'ytick',[0:0.3])
    ylabel(ax(2),'Mean Path Length.') 
    
    hold off
    
    title('Cluster Graph Clustering and Shortest Path Length')
    xlabel('Culture age')
    
    hold off;
    
    figure; hold on;
    % GF Mean cluster and shortest path
    [ax,h1,h2]=plotyy(DIVS,C_GF_mean,DIVS,dist_GF_mean,'plot');
    
    set(h1,'marker','o');
    set(h2,'marker','s');
    
    legend('C','L/N','location','best');
    
    hold(ax(1),'on');
    
    e1 = errorbar(ax(1),DIVS,C_GF_mean,errC_GF,'o');
    e1.Color = [0 0.4470 0.7410];
    set(ax(1),'YColor',[0 0 0])
%     set(ax(1),'xlim',[0 max(DIVS)],'ylim', [0,0.2],'yticklabel',[0:0.2],'ytick',[0:0.2])
    ylabel(ax(1),'Mean Clustering Coeff.')
    
    hold off;
    
    hold(ax(2),'on');
    
    e2 = errorbar(ax(2),DIVS,dist_GF_mean,errdist_GF,'s');
    e2.Color = [0.8500 0.3250 0.0980];
    set(ax(2),'YColor',[0 0 0])
%     set(ax(2),'xlim',[0 max(DIVS)],'ylim', [0,0.2],'yticklabel',[0:0.2],'ytick',[0:0.2])
    ylabel(ax(2),'Mean Path Length.') 
    
    hold off
    
    title('Full Graph Clustering and Shortest Path Length')
    xlabel('Culture age')
    
    hold off;
end