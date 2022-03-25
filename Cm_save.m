%Save data
%Create a database of the Correlation Functions
clear all
DIVS = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18];

r = exp(1:0.05:9); %Study region
for i=1:length(DIVS)
    div = DIVS(i)
    n_0 = 0;
    %Reads how many files have been recorded per DIV
    dum=sprintf('dataPLOS/*DIV%d_*.mat',DIVS(i));
    files=dir(dum);
    nf=length(files);
    Cm_DIV = zeros(5,length(r),nf); %Embbeding x length(r) x number of archives
    for n=1:nf
    %Reads the files corresponding to such DIV
        filename=horzcat(files(n).folder,'/',files(n).name);
        data=load(filename);
        
        %Net extraction
        AF=data.net.FULL_ADJACENCY; %Matriz de adyacencia completa, que incluye neuronas + bifurcaciones
        %Posición espacial de los nodos 
        xc=data.net.CLUSTER_CENTROID(:,1);
        yc=data.net.CLUSTER_CENTROID(:,2);
        %Posición epacial puntos de bifurcación
        xj=data.net.JUNCTION_CENTROID(:,1);
        yj=data.net.JUNCTION_CENTROID(:,2);
        X=[xc;xj]; 
        Y=[yc;yj];  
        
        %Buscamos el subgrafo con las componentes conexas  
        GF=graph(AF,'omitselfloops');
        %El subgrafo de la GF
%         Main cluster
%         [bin,binsize] = conncomp(GF);
%         idx = binsize(bin) == max(binsize);
%         GCF = subgraph(GF, idx);

%         Second cluster
        [bin,binsize] = conncomp(GF);
        second_component = sort(unique(binsize(bin)),'descend');
        idx = binsize(bin)  == second_component(2);
        GCF = subgraph(GF, idx);
        
        %Posicion de los nodos de la componente gigante
        X = X(find(idx));
        Y = Y(find(idx));
                    
        %Adjacency matrix for the Full Graph
        AdGCF = adjacency(GCF); %De la componente gigante
        
        %Random Walker
        n_0 = initial_node(AdGCF,n_0); %Aleatory initial node
        N = 100*numnodes(GCF); %Length of the RW
        walk = rand_walk(AdGCF,N,n_0); %Random Walker initiation
        tic
        Cm = EmbDim(5,N,walk,X,Y); %Computes the Corr. Sum 
        toc
        Cm_DIV(:,:,n) = Cm;
        

    end
    save(sprintf('Fractal_Cm_2_cluster/DIV%d_.mat',DIVS(i)'),"Cm_DIV");

end
%%
function n_0 = initial_node(AdGCF,n_0)
    if n_0 == 0 
        [row_0,~] = find(AdGCF);
        n_0 = row_0(randi(length(row_0))); %Aleatory initial node
    else
    end
    
end
function Cm = EmbDim(Emb,N,walk,X,Y)
%Algorithm to compute the Correlation Sum
    for m = 1:Emb
        j = 1;
        V = zeros(N-m,m);
        for i = 1:N-m
            V(i,:) = walk(i:i+m-1); %Serie temporal
        end
        for r=exp(1:0.05:9)
            drawnow;
            d = 0;
            heaviside = 0;
            for i = 1:N-m
                xi = X(V(i,:))';   yi = Y(V(i,:))';  
                Xj = X(V(i+1:N-m,:)); Yj = Y(V(i+1:N-m,:));
                if i == N-m-1
                    Xj = Xj'; Yj = Yj';
                end
                Xi = repmat(xi,size(Xj,1),1); Yi = repmat(yi,size(Xj,1),1);
                A = abs(Xi-Xj); B = abs(Yi-Yj);
                C = [A,B];
                aux = r-max(C');
                heaviside = heaviside+sum(aux >= 0);
            end
            Cm(m,j) = 2*heaviside/((N-m)*(N-m+1));
            j = j + 1;
        end
    end
end