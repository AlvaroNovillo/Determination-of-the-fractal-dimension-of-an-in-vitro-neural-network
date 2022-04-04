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
    allGCFdegs_2 = [];
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
        
        %Componente gigante de segundo orden 
        [bin,binsize] = conncomp(GF);
        second_component = sort(unique(binsize(bin)),'descend');
        idx_2 = binsize(bin)  == second_component(2);
        GCF_2 = subgraph(GF, idx_2);
        N_GCC_2(n,i) = numnodes(GCF_2); 
        N_GCC_2(N_GCC_2 == 0) = nan;
        
                    
        %Adjacency matrix for the Full subgraph
        AdGCF = adjacency(GCF); 
        N_GCC(n,i) = numnodes(GCF); 
        N_GCC(N_GCC == 0) = nan;
        
       
       
        %Node degree distribution
        d_GF = degree(GF);
        d_GCF = degree(GCF);
        d_GCF_2 = degree(GCF_2);
        
        %Storage of the values 
        
        [d_GF, allGFdegs] = degree_store(d_GF, allGFdegs,n);
        [d_GCF, allGCFdegs] = degree_store(d_GCF, allGCFdegs,n);
        [d_GCF_2, allGCFdegs_2] = degree_store(d_GCF_2, allGCFdegs_2,n);
        
        
        allGFdegs  = [allGFdegs d_GF];
        allGCFdegs = [allGCFdegs d_GCF];
        allGCFdegs_2 = [allGCFdegs_2 d_GCF_2];

    end
    mean_d_full(i) =mean(mean(allGFdegs,2,'omitnan'));
    mean_d_giant(i) =mean(mean(allGCFdegs,2,'omitnan'));
    mean_d_giant_2(i) =mean(mean(allGCFdegs_2,2,'omitnan'));
    

    if ismember(DIVS(i),Cases)
        div = DIVS(i)
        [Full,xf{k},pkf{k}] = degree_plot(mean(allGFdegs,2,'omitnan'), 'Node degree Full Network');
        [Giant,xg{k},pkg{k}] = degree_plot(mean(allGCFdegs,2,'omitnan'), 'Node degree GCC');
        [Giant_2,xg_2{k},pkg_2{k}] = degree_plot(mean(allGCFdegs_2,2,'omitnan'), 'Node degree GCC2');
        k = k+1;
    end  
end

%%

deg_dist(xf,pkf,'Degree distribution Full Network',Cases)
deg_dist(xg,pkg,'Degree distribution GCC',Cases)
deg_dist(xg_2,pkg_2,'Degree distribution GCC2',Cases)


%%
figure();
plot(DIVS,mean_d_full,'-o')
hold all;
plot(DIVS,mean_d_giant,'-s')
plot(DIVS,mean_d_giant_2,'g-x')
legend('Full network','GGC','GCC2','Location',"southeast")
title('Mean Degree distribution')
xlabel('DIVS')
ylabel('k')
hold off;
%%
figure();
semilogy(DIVS,mean(N_GF,1,"omitnan"),'-o')
hold on;
semilogy(DIVS,mean(N_GCC,1,"omitnan"),'-s')
semilogy(DIVS,mean(N_GCC_2,1,"omitnan"),'g-x')
legend('Full network','GCC','GCC2','Location',"southeast")
title('Node number evolution')
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
        loglog(x{1,h}(1:end-1),pk{1,h},'o')
        
    end
    hold off;
    xlabel('K');
    ylabel('P_k');
    title(comp)
    legend(strsplit(sprintf('DIV=%d ', Cases)))
    
       

end

function [d,alldegs] = degree_store(d, alldegs,n)

    if numel(d) > size(alldegs,1) && n>1
        tam=numel(d);
        alldegs(end+1:tam,:) = nan;
    else
        tam=size(alldegs,1);
        d(end+1:tam,:) = nan;
    end
end