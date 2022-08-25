function search_ann(norm_data_train,test_neurons,P,c)
% this function will test several configurations for the neural network
% and plot the surf plot correlating the number of neurons and the
% percentages of data used for training.

norm_input_data = norm_data_train(:,1:3)';
norm_target_data = norm_data_train(:,4)';

mse_train = zeros(1,1);
mse_val = zeros(1,1);    

rmse_train = zeros(1,1);
rmse_val = zeros(1,1);

fitted =  zeros(1,1);

R2_t = zeros(1,1);
R2_v = zeros(1,1);

h = 1; 
switch c
    case '1'
        for ii = test_neurons
            k = 1;
            for j = P
                %%%defining the achitecture of the ANN
                hiddenLayerSize = ii; %number of neurons - A vector defines the number of neurons per layer
                net = fitnet(hiddenLayerSize);
                net.divideParam.trainRatio = j ;
                net.divideParam.valRatio = 1-j;
                net.divideParam.testRatio = 0;
                net.trainParam.min_grad = 1e-7;
                net.trainParam.epochs = 2000;
                net.trainParam.max_fail = 6;
                
                % training the ANN
                
                [net,tr] = train( net, norm_input_data, norm_target_data);
                
                %determine the performance of the ANN
                
                modelyTrain = net(norm_input_data(:,tr.trainInd)); %%%% What is yTrain?? it's the estimated value
                modelyVal = net(norm_input_data(:,tr.valInd));
                sampleyTrain = norm_target_data(tr.trainInd); %true values
                sampleyVal = norm_target_data(tr.valInd); %true values
                
                %We don't want these values to be zero but as close to zero. 0 in training meaning most probably a sample overfitting.
                mse_train(k,h) = mean((sampleyTrain-modelyTrain).^2); % mean square error of normalized data
                mse_val(k,h) = mean((sampleyVal-modelyVal).^2); % mean square error of normalized data
                rmse_train(k,h) = sqrt(mean((sampleyTrain-modelyTrain).^2)); % mean square error of normalized data
                rmse_val(k,h) = sqrt(mean((sampleyVal-modelyVal).^2)); % mean square error of normalized data
                
                fitted (k,h) = rmse_train(k,h)./rmse_val(k,h);
                %%% as close to 1 as possible
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%% Plotter of correlation(r^2) for each configuration tested %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                ALLR2_t = corrcoef(sampleyTrain,modelyTrain);
                ALLR2_v = corrcoef(sampleyVal,modelyVal);
                
                R2_t(h)= ALLR2_t(1,2);
                R2_v(h)= ALLR2_v(1,2);
                
                R2plot(sampleyTrain,modelyTrain,sampleyVal,modelyVal,R2_t,R2_v,h);
                
                k = k+1;
            end
            h=h+1;
        end
        %grafico da relação entre dados reais e dados previstos pelo modelo
        plotterror(test_neurons,P,fitted,rmse_train,rmse_val,c);
        
    case '2'
        
        for ii = test_neurons
            k = 1;
            
            for l = 1 : P
                %%%defining the achitecture of the ANN
                hiddenLayerSize = ii; %number of neurons - A vector defines the number of neurons per layer
                net = fitnet(hiddenLayerSize);
                net.divideParam.trainRatio = 0.85 ;
                net.divideParam.valRatio = 0.15 ;
                net.divideParam.testRatio = 0;
                net.trainParam.min_grad = 1e-7;
                net.trainParam.epochs = 2000;
                net.trainParam.max_fail = 6;
                
                % training the ANN
                
                [net,tr] = train( net, norm_input_data, norm_target_data);
                
                %determine the performance of the ANN
                
                modelyTrain = net(norm_input_data(:,tr.trainInd)); %%%% What is yTrain?? it's the estimated value
                modelyVal = net(norm_input_data(:,tr.valInd));
                sampleyTrain = norm_target_data(tr.trainInd); %true values
                sampleyVal = norm_target_data(tr.valInd); %true values
                
                %We don't want these values to be zero but as close to zero. 0 in training meaning most probably a sample overfitting.
                mse_train(k,h) = mean((sampleyTrain-modelyTrain).^2); % mean square error of normalized data
                mse_val(k,h) = mean((sampleyVal-modelyVal).^2); % mean square error of normalized data
                rmse_train(k,h) = sqrt(mean((sampleyTrain-modelyTrain).^2)); % mean square error of normalized data
                rmse_val(k,h) = sqrt(mean((sampleyVal-modelyVal).^2)); % mean square error of normalized data
                
                fitted (k,h) = rmse_train(k,h)./rmse_val(k,h);
                %%% as close to 1 as possible
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%% Plotter of correlation(r^2) for each configuration tested %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                ALLR2_t = corrcoef(sampleyTrain,modelyTrain);
                ALLR2_v = corrcoef(sampleyVal,modelyVal);
                
                R2_t(h)= ALLR2_t(1,2);
                R2_v(h)= ALLR2_v(1,2);
                
                R2plot(sampleyTrain,modelyTrain,sampleyVal,modelyVal,R2_t,R2_v,h);
                
                k = k+1;
            end
            h=h+1;
        end
        %grafico da relação entre dados reais e dados previstos pelo modelo
        m = mean(rmse_train,1);
        n = mean(rmse_val,1);
        f = mean(fitted,1);
        
        plotterror(test_neurons,P,f,m,n,c);
    case '3'
        g = input('number of neurons of the 1st layer?');
        for ii = test_neurons
            k = 1;
            
            for l = 1 : P
                %%%defining the achitecture of the ANN
                hiddenLayerSize = [g,ii]; %number of neurons - A vector defines the number of neurons per layer
                net = fitnet(hiddenLayerSize);
                net.divideParam.trainRatio = 0.85 ;
                net.divideParam.valRatio = 0.15 ;
                net.divideParam.testRatio = 0;
                net.trainParam.min_grad = 1e-7;
                net.trainParam.epochs = 2000;
                net.trainParam.max_fail = 6;
                
                % training the ANN
                
                [net,tr] = train( net, norm_input_data, norm_target_data);
                
                %determine the performance of the ANN
                
                modelyTrain = net(norm_input_data(:,tr.trainInd)); %%%% What is yTrain?? it's the estimated value
                modelyVal = net(norm_input_data(:,tr.valInd));
                sampleyTrain = norm_target_data(tr.trainInd); %true values
                sampleyVal = norm_target_data(tr.valInd); %true values
                
                %We don't want these values to be zero but as close to zero. 0 in training meaning most probably a sample overfitting.
                mse_train(k,h) = mean((sampleyTrain-modelyTrain).^2); % mean square error of normalized data
                mse_val(k,h) = mean((sampleyVal-modelyVal).^2); % mean square error of normalized data
                rmse_train(k,h) = sqrt(mean((sampleyTrain-modelyTrain).^2)); % mean square error of normalized data
                rmse_val(k,h) = sqrt(mean((sampleyVal-modelyVal).^2)); % mean square error of normalized data
                
                fitted (k,h) = rmse_train(k,h)./rmse_val(k,h);
                %%% as close to 1 as possible
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%% Plotter of correlation(r^2) for each configuration tested %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                ALLR2_t = corrcoef(sampleyTrain,modelyTrain);
                ALLR2_v = corrcoef(sampleyVal,modelyVal);
                
                R2_t(h)= ALLR2_t(1,2);
                R2_v(h)= ALLR2_v(1,2);
                
                R2plot(sampleyTrain,modelyTrain,sampleyVal,modelyVal,R2_t,R2_v,h);
                
                k = k+1;
            end
            h=h+1;
        end
        %grafico da relação entre dados reais e dados previstos pelo modelo
        m = mean(rmse_train,1);
        n = mean(rmse_val,1);
        f = mean(fitted,1);
        
        plotterror(test_neurons,P,f,m,n,c);

end
end