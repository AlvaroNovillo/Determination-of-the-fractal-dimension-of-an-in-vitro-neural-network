clear all
%Scaling of the networks and fitting interval visualizer.
% All DIVS studied
DIVS = [0	1	2	3	4	5	6	7	8	9	10	11	12	13	14 16 18];
%We can represent their fractal dimension with them with fractal_paint.m

r = exp(1:0.05:9); %Study region
for i=1:length(DIVS)
    div = DIVS(i)
    %Reads how many files have been recorded per DIV
    dum=sprintf('Fractal_Cm/DIV%d_.mat',DIVS(i));
    files=load(dum);
    [~,~,nf]= size(files.Cm_DIV); %Number of files per DIV
    for n = 1:nf
        archive = n
        figure();
        for m = 1:5
            Cm = files.Cm_DIV(m,:,n);
            %Filter the data
            Cm = smoothdata(Cm,'lowess',6); 
            [r_int,int] = new_filter(r,Cm);
            hold all;
            paint(r,Cm,int,DIVS) %Paints Cm(r) vs r
%             [frac_dim(m,n),delta(m,n)] = fractalfit(r_int,Cm,int); %Computes the corr.dimension as a function of m.
        end
        legend('m=1','','m=2','','m=3','','m=4','', 'm=5','Fit interval','Location','Best');
        hold off;
    end
end

%%
function paint(r,Cm,int,DIVS)
    %Representation of Cm(r) vs r
    loglog(log(r),Cm(1,1:end),'-o','MarkerSize',4)
    loglog(log(r(int)), Cm(1,int), 'b-o','MarkerSize',4)
    xlabel('r');
    ylabel('C_m(r)');
    title('Full Graph Embedding dimension');
end


%Fit
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

function [r_int,int] = new_filter(r,Cm)
    A = diff(Cm);
%     [psor,lsor] = findpeaks(A,'SortStr','descend');
    %get the last interval where nonzero elements are
    index = (A>= max(A)/2); %Threshold
    gt=find(index~=0);
    lower = min(gt);
    upper = max(gt);
    int = lower:upper; %Interval where the fit is performed
    r_int = log(r(int)); %Differential section of r where the fit is performed 
%     Visualization of the alorithim
%     figure();
%     plot(A)
%     hold on;
%     plot(int,A(int),'r-o')
%     legend('Cm differences','Fit interval')
%     hold off;
end