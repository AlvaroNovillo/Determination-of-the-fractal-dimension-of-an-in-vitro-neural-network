%Save data
%Create a database of the Correlation Functions
clear all
DIVS=[0     1     2     4     7     9    11    16    18];
n_0 = 0;
r = exp(1:0.05:9); %Study region
for i=1:length(DIVS)
    %Reads how many files have been recorded per DIV
    dum=sprintf('dataPLOS/*DIV%d_*.mat',DIVS(i));
    files=dir(dum);
    nf=length(files);
    Cm_DIV = zeros(length(r),nf+1);
    for n=1:nf
    %Reads the files corresponding to such DIV
        filename=horzcat(files(n).folder,'/',files(n).name);
        data=load(filename);
        
        %extraes la red
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
        [bin,binsize] = conncomp(GF);
        idx = binsize(bin) == max(binsize);
        GCF = subgraph(GF, idx);
        %Posicion de los nodos de la componente gigante
        X = X(find(idx));
        Y = Y(find(idx));
                    
        %Adjacency matrix for the Full Graph
        AdGCF = adjacency(GCF); %De la componente gigante
        
        %Random Walker
        rng(123); %Seed of the random walker fixed
        n_0 = initial_node(AdGCF,n_0); %Aleatory initial node
        N = numnodes(GCF); %Length of the RW
        walk = rand_walk(AdGCF,N,n_0); %Random Walker initiation
        tic
        Cm = EmbDim(5,N,walk,X,Y,r); %Computes the Corr. Sum 
        toc
        Cm_DIV(:,1) = r;
        Cm_DIV(:,n+1) = Cm';

    end
    save(sprintf('Fractal_Cm/DIV%d_.mat',DIVS(i)'),"Cm_DIV");

end
%%
function n_0 = initial_node(AdGCF,n_0)
    if n_0 == 0 
        [row_0,~] = find(AdGCF);
        n_0 = row_0(randi(length(row_0))); %Aleatory initial node
    else
    end
    
end
function Cm = EmbDim(Emb,N,walk,X,Y,r)
%Algorithm to compute the Correlation Sum
    for m = Emb
        j = 1;
        V = zeros(N-m,m);
        for i = 1:N-m
            V(i,:) = walk(i:i+m-1); %Serie temporal
        end
        for r=r
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
            Cm(1,j) = 2*heaviside/((N-m)*(N-m+1));
            j = j + 1;
        end
    end
end