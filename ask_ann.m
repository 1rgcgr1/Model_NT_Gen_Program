function  [output]=ask_ann(input,net)
%% input =[Wm,VDm,x];
    
    [norm_input_data,train_target_min,train_target_range] = norma(input,'obtain');
     
        outputs = net(norm_input_data'); 
        
        
        % denormalizing data obtained from network
        target_min = -1;
        target_range = 2;
        % denormalization expression
        for g = 1: size(outputs,2)
            den_net_out =exp((((outputs(g) - target_min)./target_range).*train_target_range)+ train_target_min);
            if g == 1
                den_out = den_net_out;
            else
                den_out = [den_out;den_net_out];
            end
        end
        output = den_out;
end