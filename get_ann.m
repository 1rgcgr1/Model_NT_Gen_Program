function [rmse_train,rmse_val,net] = get_ann(norm_data_train,test_neurons,P)

norm_input_data = norm_data_train(:,1:size(norm_data_train,2)-1)';
norm_target_data = norm_data_train(:,size(norm_data_train,2))';



%%%defining the achitecture of the ANN
hiddenLayerSize = test_neurons; %number of neurons - A vector defines the number of neurons per layer
net = fitnet(hiddenLayerSize);
net.trainFcn = 'trainlm';
net.divideParam.trainRatio = P ;
net.divideParam.valRatio = 1-P;
net.divideParam.testRatio = 0;
net.trainParam.min_grad = 1e-7;
net.trainParam.epochs = 2000;
net.trainParam.max_fail = 6;

% training the ANN

[net,tr] = train( net, norm_input_data, norm_target_data);

%determine the performance of the ANN
                                  
modelyTrain = net(norm_input_data(:,tr.trainInd)); %%%%
modelyVal = net(norm_input_data(:,tr.valInd));
sampleyTrain = norm_target_data(tr.trainInd); %true values
sampleyVal = norm_target_data(tr.valInd); %true values

%We don't want these values to be zero but as close to zero. 0 in training meaning most probably a sample overfitting.
%mse_train = mean((sampleyTrain-modelyTrain).^2); % mean square error of normalized data
%mse_val = mean((sampleyVal-modelyVal).^2); % mean square error of normalized data
rmse_train = sqrt(mean((sampleyTrain-modelyTrain).^2)); % mean square error of normalized data
rmse_val = sqrt(mean((sampleyVal-modelyVal).^2)); % mean square error of normalized data


      
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Plotter of correlation(r^2) for each configuration tested %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
       

ALLR2_t = corrcoef(sampleyTrain,modelyTrain);
ALLR2_v = corrcoef(sampleyVal,modelyVal);

R2_t = ALLR2_t(1,2);
R2_v = ALLR2_v(1,2);
%function R2plot

R2plot(sampleyTrain,modelyTrain,sampleyVal,modelyVal,R2_t,R2_v,1);


end