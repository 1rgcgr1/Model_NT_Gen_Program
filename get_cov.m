function [store_cov] = get_cov(CV_data,FixedL)

%store_cov = [Tlength,Twidth,VG_SD,Cap];
%Get CV data on a specific L
fL = FixedL;

result_L = CV_data(:,1) == FixedL*1e-6;

CV_Data_W = CV_data(result_L,2:4);

i=1;
while any(CV_Data_W)
    
    result_W = CV_Data_W(:,1) == CV_Data_W(1,1);
    
    Worg = CV_Data_W(result_W,:); % um w especifico
    ppl = (size(Worg,1)/sum( Worg(:,2) == Worg(1,2)));
    %        plotVD = train(1:(size(train(:,3),1)/sum( train(:,3)==train(1,3) )),3);
    for k = 0:(size(Worg,1)/ppl)-1
        
        %         VG = Worg(1:(size(Worg(:,2),1)/sum( Worg(:,2)==Worg(1,2) )),2); 
        if k == 0
            curves = Worg(1:ppl,3);
        else
            curves = [curves, Worg(k*ppl+1:(k+1)*ppl,3)];
        end
        
        
    end
    
    Cap = mean(curves,2); % guardado lateralmente
    fW = Worg(1,1)*1e6;
    Cov = round(min(Cap)/2*1e12,3); %in pF
    min_c = min(Cap);
    max_c = max(Cap);
    Cox = (max_c-min_c)/(fW*1e-4*fL*1e-4)*1e9 ;% in nF/cm^2;
    
    % Store information
    if i == 1
        store_cov = [fW,Cov,Cox];
    else
        store_cov = [store_cov;fW,Cov,Cox];
    end
    
    %Get individual values of Cov and Cox for each Ws
    i = i+1;
    CV_Data_W = CV_Data_W(~result_W,:);
end
end