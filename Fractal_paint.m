clear all
%Study of the evolution of the fractal dimension per DIV
DIVS = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18];

for h= 1:2 %Clusters loop: loops in giant components
    r = exp(1:0.05:9); %Study region
    for i=1:length(DIVS)
        %Reads how many files have been recorded per DIV
        dum=sprintf('Fractal_Cm_%d/DIV%d_.mat',h,DIVS(i));
        files=load(dum);
        [~,~,nf]= size(files.Cm_DIV); %Number of files per DIV
        for n = 1:nf
            for m = 1:5
                Cm = files.Cm_DIV(m,:,n);
                %Filter the data
                Cm = smoothdata(Cm,'lowess',6); 
                [r_int,int] = new_filter(r,Cm);
                [frac_dim(m,n),delta(m,n)] = fractalfit(r_int,Cm,int); %Computes the corr.dimension as a function of m.
            end
        end
        %Mean of the fractal dimension per m
        if nf > 1
            frac_dim(:,nf+1) = mean(frac_dim,2);
            delta(:,nf+1) = sqrt(mean(delta,2).^2 + (tinv(.95,nf)*sqrt(1/(nf-1)*sum(sum(delta-mean(delta,2))).^2)).^2);
            mean_frac_dim(:,i) =  frac_dim(:,nf+1);
            mean_delta(:,i) = delta(:,nf+1);
        else
            mean_frac_dim(:,i) =  frac_dim(:,nf);
            mean_delta(:,i) = delta(:,nf);
        end
    end
    figure();
    hold all;
    for k = 1:5
        frac_paint(mean_frac_dim(k,:),mean_delta(k,:), DIVS)
    end
    legend('m=1','m=2','m=3','m=4', 'm=5','Location','Best');
    title(sprintf('Giant component %d'),h)
    hold off;
    frac_m5(h,:) = mean_frac_dim(5,:);
    delta_frac_m5(h,:) = mean_delta(5,:);
end
%%
figure();
hold on;
arrayfun(@(x) errorbar(DIVS,frac_m5(x,:),delta_frac_m5(x,:),'-o'),1:2)
xlabel('DIVs')
ylabel('\beta') 
legend('GCC','GCC2')
%%
function [frac_dim,delta] = fractalfit(r_int,Cm,int)
    %Computation of the fit
    [P_5,S] = polyfit(r_int,log(Cm(1,int)),1);
    frac_dim = P_5(1);
    %Estimation of the standard error of the slope
    delta = sqrt(1/(length(r_int)-2)*(sum((log(Cm(1,int)) - mean(log(Cm(1,int)))).^2))/sum((r_int - mean(r_int)).^2));
    %Linear fit statistics
    Rscore = 1 - (S.normr/norm(log(Cm(1,int)) - mean(log(Cm(1,int)))))^2;


end

function [r_int,int] = new_filter(r,Cm)
    A = diff(Cm);
    index = (A>= max(A)/2); %Threshold
    gt=find(index~=0);
    lower = min(gt);
    upper = max(gt);
    int = lower:upper; %Interval where the fit is performed
    r_int = log(r(int)); %Differential section of r where the fit is performed 
end

function frac_paint(frac_dim,delta, DIVS)
    errorbar(DIVS, frac_dim,delta, '-o')
%     legend('Giant component corr. dimension')
    xlabel('DIVs')
    ylabel('\beta')   
end