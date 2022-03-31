%Node degree distribution D = degree(g)
clear all
DIVS = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18];
k = 1;
Cases = input('Enter cases with [ ] around them');
for i=1:length(DIVS)
    dum=sprintf('dataPLOS/*DIV%d_*.mat',DIVS(i));
    files=dir(dum);
    nf=length(files);
    allGFdegs = [];
    allGCFdegs = [];
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
        N_GF(n,i) = numnodes(GF); 
        N_GF(N_GF == 0) = nan;
        
        %El subgrafo de la GF
        [bin,binsize] = conncomp(GF);
        idx = binsize(bin) == max(binsize);
        GCF = subgraph(GF, idx);
                    
        %Adjacency matrix for the Full subgraph
        AdGCF = adjacency(GCF); 
        N_GCC(n,i) = numnodes(GCF); 
        N_GCC(N_GCC == 0) = nan;
       
       
        %Node degree distribution
        d_GF = degree(GF);
        d_GCF = degree(GCF);
        
        %Storage of the values 
        
        if numel(d_GF) > size(allGFdegs,1) && n>1
            tam=numel(d_GF);
            allGFdegs(end+1:tam,:) = nan;
        else
            tam=size(allGFdegs,1);
            d_GF(end+1:tam,:) = nan;
        end

        if numel(d_GCF) > size(allGCFdegs,1) && n>1
            tam=numel(d_GCF);
            allGCFdegs(end+1:tam,:) = nan;
        else
            tam=size(allGCFdegs,1);
            d_GCF(end+1:tam,:) = nan;  
        end
        allGFdegs  = [allGFdegs d_GF];
        allGCFdegs = [allGCFdegs d_GCF];

    end
    mean_d_full(i) =mean(mean(allGFdegs,2,'omitnan'));
    mean_d_giant(i) =mean(mean(allGCFdegs,2,'omitnan'));
    

    if ismember(DIVS(i),Cases)
        div = DIVS(i)
        [Full,xf{k},pkf{k}] = degree_plot(mean(allGFdegs,2,'omitnan'), 'Node degree GFull');
        [Giant,xg{k},pkg{k}] = degree_plot(mean(allGCFdegs,2,'omitnan'), 'Node degree GGiant');
        k = k+1;
    end  
end

%%

deg_dist(xf,pkf,'Degree distribution GFull',Cases)
deg_dist(xg,pkg,'Degree distribution GGiant',Cases)


%%
figure();
plot(DIVS,mean_d_full,'-o')
hold all;
plot(DIVS,mean_d_giant,'-s')
legend('Full network','GGC')
title('Mean Degree distribution')
xlabel('DIVS')
ylabel('k')
hold off;
%%
figure();
semilogy(DIVS,mean(N_GF,1,"omitnan"),'-o')
hold on;
semilogy(DIVS,mean(N_GCC,1,"omitnan"),'-s')
legend('Full network','GGC')
title('Nodes GGC')
xlabel('DIVS')
ylabel('Nodes')
hold off
%%
function [plot,x,pk] = degree_plot(D,comp)
    figure("Visible","off");
    h = histogram(D,unique(D));
    x = h.BinEdges;
    y = h.Values;
    pk = y/length(D);
    
    figure();
    plot = semilogy(x(1:end-1), y,'o');
    xlabel('Node degree');
    ylabel('Node number');
    title(comp) 
end
function deg_dist(x,pk,comp,Cases)
    figure();
    
    hold all;
    for h = 1:length(x) 
        plot(x{1,h}(1:end-1),pk{1,h},'o')
        
    end
    hold off;
    xlabel('K');
    ylabel('P_k');
    title(comp)
    legend(strsplit(sprintf('DIV=%d ', Cases)))
    
       

end