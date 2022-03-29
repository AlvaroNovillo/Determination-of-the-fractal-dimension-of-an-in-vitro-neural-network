%Node degree distribution D = degree(g)
clear all
DIVS = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18];

%r= rango scaling
for i=1:length(DIVS)
    dum=sprintf('dataPLOS/*DIV%d_*.mat',DIVS(i));
    files=dir(dum);
    nf=length(files);

    for n=1:nf

        filename=horzcat(files(n).folder,'/',files(n).name);
        data=open(filename);
        %extraes la red
        AF=data.net.FULL_ADJACENCY; %Matriz de adyacencia completa, que incluye neuronas + bifurcaciones
        
        %Posición espacial de los nodos 
        xc=data.net.CLUSTER_CENTROID(:,1);
        yc=data.net.CLUSTER_CENTROID(:,2);
        
        %Buscamos el subgrafo con las componentes conexas  
        GF=graph(AF,'omitselfloops');
        %El subgrafo de la GF
        [bin,binsize] = conncomp(GF);
        idx = binsize(bin) == max(binsize);
        GCF = subgraph(GF, idx);
                    
        %Adjacency matrix for the Full subgraph
        AdGCF = adjacency(GCF); 
        N(n,i) = numnodes(GCF); 
        N(N == 0) = nan;
        %Preparation of the arrays
        [row_G,~] = size(AdGCF); %Rows
        [row_F,~] = size(AF);

       
        %Node degree distribution
        D_Full(1:row_F,n) = degree(GF);
        D_Giant(1:row_G,n) = degree(GCF);
    end
    div = DIVS(i)
    mean_d_full(i) =mean(sum(D_Full,2))/length(sum(D_Full,2));
    mean_d_giant(i) =mean(sum(D_Giant,2))/length(sum(D_Giant,2));
    Full = degree_plot(sum(D_Full,2), 'Node degree Full');
    Giant = degree_plot(sum(D_Giant,2), 'Node degree Giant');
end
%%
figure();
semilogy(DIVS,mean_d_full,'-o')
hold all;
semilogy(DIVS,mean_d_giant,'-s')
legend('Full network','GGC')
title('Mean Degree distribution')
xlabel('DIVS')
ylabel('k')
hold off
%%
figure();
semilogy(DIVS,mean(N,1,"omitnan"),'-o')
title('Nodes GGC')
xlabel('DIVS')
ylabel('Nodes')
%%
function plot = degree_plot(D,comp)
    figure("Visible","off");
    h = histogram(D,unique(D));
    x = h.BinEdges;
    y = h.Values;
    
    figure();
    plot = semilogy(x(1:end-1), y,'o');
    xlabel('Node degree');
    ylabel('Node number');
    title(comp)
    
    
    
    
    
end