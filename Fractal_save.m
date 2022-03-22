clear all
%Analysis of the data

DIVS = [0	1	2	3	4	5	6	7	8	9	10	11	12	13	14 16 18];
%Parameters obtained from the analysis
%We can represent them with fractal_paint.m
% k = [4 2 4 3 4 5 6 4 6 4 6 4 9 4 6 7 2]; %Order of the local maxima
% n = [4 3 3 3 7 2 3 6 2 5 2 2 2 2 2 2 2]; %File we want to take

DIVS = 5;
for i=1:length(DIVS)
    %Reads how many files have been recorded per DIV
    dum=sprintf('Fractal_Cm/DIV%d_.mat',DIVS(i));
    files=load(dum);
    r = files.Cm_DIV(:,1); %Intervalo de estudio
    [~,nf]= size(files.Cm_DIV); %numero de archivos
    for n = 2:nf
        n
        Cm = files.Cm_DIV(:,n)';
        %Filter the data
        Cm = smoothdata(Cm,'lowess',6); 
        [r_int,int] = new_filter(r,Cm,4);
        paint(r,Cm,int) %Paints Cm(r) vs r
        frac_dim = fractalfit(r_int,Cm,int) %Computes the corr.dimension as a function of m.
        
    end
end
%%
function paint(r,Cm,int)
    %Representation of Cm(r) vs r
    figure();
    loglog(r,Cm(1,1:end),'g-o','MarkerSize',4)
    hold all;
    loglog(r(int), Cm(1,int), 'b-o','MarkerSize',4)
    xlabel('r');
    ylabel('C_m(r)');
    title('Full Graph Embedding dimension');
    legend('m=5','fit interval',"Location",'Best');
    hold off;
    
end


function slope = PartialSlopes(r,Cm,h)
    %Compute several linear fits in our graphs    
    % %Arguments: r ->Study region
    %           Cm ->Correlation Sum
    %           h-> Array denoting the intervals where the fits are going
    %           to be made
    
    slope = zeros(1,length(h));
    for g = 1:length(h)-1
        r_int = log(r(h(g):h(g+1)));
        P_5 = polyfit(r_int,log(Cm(1,h(g):h(g+1))),1);
        slope(g) = P_5(1);
    end
    slope(slope>2.5) = 0; %Filter values of the slope > 2
    slope(slope<0) = 0; %Filter negative slope values
    
end


%Fit
function frac_dim = fractalfit(r_int,Cm,int)
    %Computation of the fits 
%     P_1 = polyfit(r_int,log(Cm(1,int)),1);
%     P_2 = polyfit(r_int,log(Cm(2,int)),1);
%     P_3 = polyfit(r_int,log(Cm(3,int)),1);
%     P_4 = polyfit(r_int,log(Cm(4,int)),1);
    [P_5,S] = polyfit(r_int,log(Cm(1,int)),1);
%     beta = [P_1(1),P_2(1),P_3(1),P_4(1),P_5(1)];
    %Plotting \beta against r, The correlation exponent increases
    %linearly with the embedding dimension m
%     hold on;
%     figure();
%     plot(1:5,beta,'o-')
%     xlabel('m');
%     ylabel('\beta'); 
%     legend('\beta');
%     hold off;
    frac_dim = P_5(1);
    Rscore = 1 - (S.normr/norm(log(Cm(1,int)) - mean(log(Cm(1,int)))))^2;
    sprintf('The Correlation Dimension is %.3f, with R^2 of %.3f',P_5(1),Rscore)
    


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