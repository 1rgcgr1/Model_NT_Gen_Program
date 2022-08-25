function [ norm_data,target_min,target_range] = norma(data_train,c)
%Expression used: (target_range*((input_data(:,i)-train_input_min(i))/train_input_range(i)))+target_min; %function minmax
        
%preprocess data so that minimum is -1 and maximum is 1
        
% for reference [W,VD,VG]
norma_range =  2 ; %%to keep the data values in [-1 1] range
norma_min = -1;

switch c
    case 'define'
        %Processing data_train
        %Defining values for max, min, range.
        
        train_min = min(data_train);
        train_max = max(data_train);
        train_range = train_max - train_min;
        
        target_min = train_min(1,size(train_min,2));
        target_range = train_range(1,size(train_range,2));
        norma_train = [train_min;train_max;train_range] ;
        
        save(fullfile(tempdir, 'norma_data.mat'), 'norma_train', '-mat');
        
    case 'obtain'
        FileData = load(fullfile(tempdir, 'norma_data.mat'));
        norma_Data = FileData.norma_train;
        
        train_min = norma_Data(1,:);
        train_max = norma_Data(2,:);
        train_range = norma_Data(3,:);
        
        target_min = train_min(1,size(train_min,2));
        target_range = train_range(1,size(train_range,2));
       
        %use previously defined values for max and min used to define the
        %model. In order to obtain the normalized values of a test matrix
        %it can change the input values but it still has to take
        %into account the min an max previously obtained previously.
end
norm_data = zeros(size(data_train))';


for i =1:size(data_train,2)
    aux_norm_data =(norma_range*((data_train(:,i)-train_min(i))/train_range(i)))+norma_min; %function minmax between -1 and 1
    if i == 1
        norm_data = aux_norm_data ;
    else
        norm_data = [norm_data,aux_norm_data];
    end
    
end



%%%%%% to avoid errors uncomment 
norm_data(isnan(norm_data))= 0;
% norm_target_data(isnan(norm_target_data))=0;

    
end