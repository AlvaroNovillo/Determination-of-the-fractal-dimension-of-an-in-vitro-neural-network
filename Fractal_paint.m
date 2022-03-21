clear all
%Results
DIVS=[0     1     2     4     7     9    11    16    18];
%Parameters for the filter
k = [4 2 4 4 4 4 4 7 2];%Order of the local maxima
n = [4 3 2 7 6 5 2 2 2]; %File we want to take

for i=1:length(DIVS)
    DIVS(i)
    %Reads how many files have been recorded per DIV
    dum=sprintf('Fractal_Cm/DIV%d_.mat',DIVS(i));
    files=load(dum);
    r = files.Cm_DIV(:,1); %Intervalo de estudio
    [~,nf]= size(files.Cm_DIV); %numero de archivos
    Cm = files.Cm_DIV(:,n(i))';
    %Filter the data
    Cm = smoothdata(Cm,'lowess',6); 
    [r_int,int] = new_filter(r,Cm,k(i));
    %Results
    paint(r,Cm,int) %Paints Cm(r) vs r
    [frac_dim(i),delta(i)] = fractalfit(r_int,Cm,int); %Returns the corr.dim per DIV.

end
frac_paint(frac_dim,delta,DIVS) %Returns the graph of the evolution of beta
%%
function paint(r,Cm,int)
    %Representation of Cm(r) vs r
    figure();
    loglog(r,Cm(1,1:end),'-o','MarkerSize',4)
    hold all;
    loglog(r(int), Cm(1,int), 'r-o','MarkerSize',4)
    xlabel('r');
    ylabel('C_m(r)');
    title('Full Graph Embedding dimension');
    legend('m=5','fit interval',"Location",'Best');
    hold off;
    
end




function [frac_dim,delta] = fractalfit(r_int,Cm,int)
    %Computation of the fit
    [P_5,S] = polyfit(r_int,log(Cm(1,int)),1);
    frac_dim = P_5(1);
    %Estimation of the standard error of the slope
    delta = sqrt(1/(length(r_int)-2)*(sum((log(Cm(1,int)) - mean(log(Cm(1,int)))).^2))/sum((r_int - mean(r_int)).^2));
    %Linear fit statistics
    Rscore = 1 - (S.normr/norm(log(Cm(1,int)) - mean(log(Cm(1,int)))))^2;
    sprintf('The Correlation Dimension is %.3f +/- %.3f, with R^2 of %.3f',P_5(1),delta,Rscore)
    
end

function [r_int,int] = new_filter(r,Cm,k)
    A = diff(Cm);
    [psor,lsor] = findpeaks(A,'SortStr','descend');
    %get the last interval where nonzero elements are
    index = (A>= psor(k)); %Only greater than the kth local max
    gt=find(index~=0);
    lower = min(gt)-1;
    upper = max(gt)-2;
    int = lower:upper; %Interval where the fit is performed
    r_int = log(r(int)); %Differential section of r where the fit is performed 
    figure();
    hold on;
    findpeaks(A)
    plot(int,A(int),'r-o')
    legend('Cm differences', 'Local maxima','Fit interval')
    hold off;
end

function frac_paint(frac_dim,delta, DIVS)
    figure();
    errorbar(DIVS, frac_dim,delta, '-o')
    legend('Giant component corr. dimension')
    xlabel('DIVs')
    ylabel('\beta')
    
    
end