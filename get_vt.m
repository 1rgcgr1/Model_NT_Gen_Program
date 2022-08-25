function [store_VT] = get_vt(Data,FixedL,c)

%store_vt = [Tlength,Twidth,VG,ID];
test_neurons = 5;
P = 0.8;

fL = FixedL;

result_L = Data(:,1) == FixedL*1e-6;

Data_W = Data(result_L,2:end); %% [W,VG,ID];

i=1;
while any(Data_W)
    
    result_W = Data_W(:,1) == Data_W(1,1);
    
    Worg = Data_W(result_W,:); % um w especifico
    ppl = (size(Worg,1)/sum( Worg(:,2) == Worg(1,2)));
%            plotVD = train(1:(size(train(:,3),1)/sum( train(:,3)==train(1,3) )),3);
    for k = 0:(size(Worg,1)/ppl)-1 %Separating the similar curves
        
        %         VG = Worg(1:(size(Worg(:,2),1)/sum( Worg(:,2)==Worg(1,2) )),2);
        if k == 0
            curves = Worg(1:ppl,3);
            VG = Worg(1:ppl,2);
        else
            curves = [curves,Worg(k*ppl+1:(k+1)*ppl,3)];
        end
        
    end
    switch c
        case 'lin'
            
            curv = mean(curves,2); 
        case 'sat'
            curv = sqrt(mean(curves,2));
    end
    fW = Worg(1,1)*1e6;
    
    data_train = [VG,curv];
    
    [norm_data_train,train_target_min,train_target_range] = norma(data_train,'define');
    [rmse_train,rmse_val,net] = get_ann(norm_data_train,test_neurons,P);
%         plot_model_data_BB(net,plotVG,IDcurve,norm_input_data,norm_target_data,train_target_min,train_target_range,fW_l,fL_l);
    
    outputs = net(norm_data_train(:,1)');
    
    % denormalizing data obtained from network
    target_min = -1;
    target_range = 2;
    % denormalization expression
    for g = 1: size(outputs,2)
        den_net_out =exp(((outputs(g) - target_min)./target_range).*train_target_range)+ train_target_min;
       if g == 1
           den_out = den_net_out;
       else
           den_out = [den_out;den_net_out];
       end
    end
    
    nd_gm_m = diff(diff(den_out));
    
    ind = find(nd_gm_m == max(nd_gm_m),1);
    VT = VG(ind);
    
    % Store information
    if i == 1
        store_VT = [fW,fL,VT];
    else
        store_VT = [store_VT;fW,fL,VT];
    end
    
    %Get individual values of Cov and Cox for each Ws
    i = i+1;
    Data_W = Data_W(~result_W,:);
end
end


