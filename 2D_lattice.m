clear all
%Creation of (10,10) 2D lattice
G=net_lattice(10);
A = adjacency(G);
m = ones(10,10);
[X,Y] = find(m);
plot(G,'XData',X,'YData',Y)
%In order to see clearly the 2D lattice that we are representing
%%
n_0 = 0;
r = -5:15;
%Random Walker
n_0 = initial_node(A,n_0); %Aleatory initial node
N = 2000; %Length of the RW
walk = rand_walk(A,N,n_0); %Random Walker initiation
tic
Cm = EmbDim(4,N,walk,X,Y); %Computes the Corr. Sum 
toc
%%
%Plot Cm(r) vs r
[r_int,int,A] = arrayfun(@(x) new_filter(r,Cm(x,:)),1:4,'UniformOutput',false);
figure();
hold on;
arrayfun(@(x) paint(r,Cm(x,:),int{x}),1:4)
legend('m=1','','m=2','','m=3','','m=4','Fit interval','Location','Best');
hold off;
%Visualization of the algorithim
figure();
plot(A{4})
hold on;
plot(int{4},A{4}(int{4}),'r-o')
legend('Cm differences (m=4)','Fit interval')
hold off;
%Computation of the fractal dimension 
[frac_dim,~] = arrayfun(@(x) fractalfit(r_int{x},Cm(x,:),int{x}),1:4,'UniformOutput',false)
%\beta vs m
figure();
hold on;
arrayfun(@(x) plot(x,frac_dim{x},'ko'), 1:4)
hold off
xlabel('m');
ylabel('\beta'); 
legend('\beta');
hold off;

%%
function paint(r,Cm,int)
    %Representation of Cm(r) vs r
    loglog(r,Cm(1,1:end),'-o','MarkerSize',4)
    loglog(r(int), Cm(1,int), 'b-o','MarkerSize',4)
    xlabel('r');
    ylabel('C_m(r)');
    title('Full Graph Embedding dimension');
end
function [r_int,int,A] = new_filter(r,Cm)
    A = diff(Cm);
    index = (A>= max(A)/2); %Threshold
    gt=find(index~=0);
    lower = min(gt);
    upper = max(gt);
    int = lower:upper; %Interval where the fit is performed
    r_int = log(r(int)); %Differential section of r where the fit is performed 

end

function [frac_dim,delta] = fractalfit(r_int,Cm,int)
    %Computation of the fit
    [P_5,S] = polyfit(r_int,log(Cm(1,int)),1);
    frac_dim = P_5(1);
    %Estimation of the standard error of the slope
    delta = sqrt(1/(length(r_int)-2)*(sum((log(Cm(1,int)) - mean(log(Cm(1,int)))).^2))/sum((r_int - mean(r_int)).^2));
    %Linear fit statistics
    Rscore = 1 - (S.normr/norm(log(Cm(1,int)) - mean(log(Cm(1,int)))))^2;
%     sprintf('The Correlation Dimension is %.3f +/- %.3f, with R^2 of %.3f',P_5(1),delta,Rscore)


end

function Cm = EmbDim(Emb,N,walk,X,Y)
%Algorithm to compute the Correlation Sum
    for m = 1:Emb
        j = 1;
        V = zeros(N-m,m);
        for i = 1:N-m
            V(i,:) = walk(i:i+m-1);
        end
        for r = -5:15
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


function n_0 = initial_node(AdGCF,n_0)
    if n_0 == 0
        [row_0,~] = find(AdGCF);
        n_0 = row_0(randi(length(row_0))); %Aleatory initial node
    else
    end
end

function [G]=net_lattice(n)
A=delsq(numgrid('S',n+2));
G=graph(A,'omitselfloops');
end


function walk = rand_walk(R,N,n)
% R --> matriz de adyacencia
% N --> numero de elementos del random walker
% n --> nodo inicial
    walk = zeros(1,N);
    for i = 1:N
        neigh = find (R(:,n));
        walk(i) = neigh(randi(length(neigh)));
        n = walk(i);
    end

end