% 1.1st version
% Author: Rui Guilhoto
% Intrusctions: 
% Warning: this program does not contain the required code to run on all computers 
% That functionality can be added later. For now the program can run entirely on command window on my pc
% as long as data is structured in similar fashion.
%
% In order to use the menus, press the number you would like to use and then
% press enter. In order for the program to properly work, the information
% must be colected in the lab using always the same program and the same structure in label. 
% Any change on the program may require additional configuration 
%
% The program is stuctured in a linear fashion but allows users to jump
% functionalities.
% [1] - This step is always required to run first. It gets the raw data
% from the original lab files and stores a temporary file containing the
% required information for the other steps to work. It also allows users to
% check for abnormal data(outliers) and check for variability from the
% collected samples. In this step you should be looking for the most
% reliable aspect ratios. this step will also store the generated figures.
%
% [2] - This step allows users to test a few configurations on the ANN
% model and check if it works properly for the previously specified L.
%
% [3] - This step usually occurs after step 2 and will generate the ANN model
%
% [4] - This step will allow to verify if the new generated model is predicting well the behavior of the transistors
% 
% [5] - This step will allow the user to test the first derivative of the model to check for model continuity 
%
% [6] - This step will generate the verilogA file required to run in Cadence
%
% [7] - Close the program
%
% For more information
% please contact me by email: r.guilhoto@campus.fct.unl.pt

function Model_program1_1
clc;

clear res
clear prompt_s
clear prompt_s1

close all;
format long
cd('D:\universidade\5ano\2semestre_tese\MAtlab'); % workspace


%go to workspace and select add to path all folders and subfolders

res = '';
prompt_s = input('&&&&&&&&&&& MENU &&&&&&&&&&&&&&&&\n[1] - Obtain data/ Manual outlier detection/ Statistical information\n[2] - Test different numbers of neurons \n[3] - Create ANN model \n[4] - Plot Model vs transistors caratheristics curves \n[5] - Check model vs transistors transconductance(gm)and output conductance(gds)\n[6] - Save model and generate verilog A file\n[7] - Close Program \n Which program do you want to run [1,2,3,4,5,6,7]?','s');
clc
switch prompt_s
    case '1' %%%%% Obtain data/ Manual outlier detection/ Statictical information
        
        for  lol = 1
            disp('&&&&&&&&&&& Sub-MENU 1 &&&&&&&&&&&&&&&&');
            disp('What data would you like to check');
            prompt_s1 = input('[1] - Saturation\n[2] - Linear  \n[3] - Output\n[4] - CV\n[5] - Define L\n[6] - ALL\n[7] - Close program\n Which one would you like to run [1,2,3,4,5,6,7]?','s');
            clc
            disp('What data would you like to check');
            switch prompt_s1
                case '1'
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get saturation IV curves %%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    clc
                    
                    [sat_train,VD_sat] = get_sat_data1_1;
                    
                    %        sat_train = [ Tlength ,Twidth, VG, MeanAbsID, stdAbsID , number of samples ]; %current is in uA
                    save(fullfile(tempdir, 'VD_sat.mat'), 'VD_sat', '-mat');
                   
                    save(fullfile(tempdir, 'sat_data.mat'), 'sat_train', '-mat');
                    
                case '2'
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get linear IV curves %%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    clc
                    
                    [lin_train,VD_lin] = get_lin_data1_1;
                    
                    %        lin_train = [ Tlength ,Twidth, VG, MeanAbsID, stdAbsID , number of samples  ]; %current is in uA
                    save(fullfile(tempdir, 'VD_lin.mat'), 'VD_lin', '-mat');
                    
                    save(fullfile(tempdir, 'lin_data.mat'), 'lin_train', '-mat');
                    
                case '3'
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get output IV curves %%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    clc
                    
                    [out_train] = get_out_data1_2 ;
                    
                    %       out_train = [ Tlength ,Twidth, VD, VG, MeanAbsID, stdAbsID , number of samples ]; %current is in uA
                    
                    save(fullfile(tempdir, 'output_data.mat'), 'out_train', '-mat');
                    
                case '4'
                    clc
                    
                    [tr_CV_train] = get_tCV_data1_1 ;
                    
                    %         tr_CV_train = [ Tlength , Twidth, V, MeanC, stdC, number of samples];
                  
                    save(fullfile(tempdir, 'CV_data.mat'), 'tr_CV_train', '-mat');
                case '5'
                    clc
                    
                    Fixed_L = input(['What is the transistor length you would like to run the program on[', char(181), 'm] ? ']);
                    save(fullfile(tempdir, 'Fixed_L.mat'), 'Fixed_L', '-mat');
                    
                case '6'
                    clc
                    
                    [sat_train,VD_sat] = get_sat_data1_1;
                    
                    save(fullfile(tempdir, 'VD_sat.mat'), 'VD_sat', '-mat');
                    save(fullfile(tempdir, 'sat_data.mat'), 'sat_train', '-mat');
                    
                    clc
                    
                    [lin_train,VD_lin] = get_lin_data1_1;
                    
                    save(fullfile(tempdir, 'VD_lin.mat'), 'VD_lin', '-mat');
                    save(fullfile(tempdir, 'lin_data.mat'), 'lin_train', '-mat');
                    
                    
                    clc
                    
                    [out_train] = get_out_data1_2 ;
                    
                    save(fullfile(tempdir, 'output_data.mat'), 'out_train', '-mat');
                    
                    
                    clc
                    
                    [tr_CV_train] = get_tCV_data1_1 ;
                   
                    save(fullfile(tempdir, 'CV_data.mat'), 'tr_CV_train', '-mat');
                    
                    clc
                    
                    Fixed_L = input(['What is the transistor length you would like to run the program on[', char(181), 'm] ? ']);
                    save(fullfile(tempdir, 'Fixed_L.mat'), 'Fixed_L', '-mat');
                    
                case '7'
                    clc
                    return
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            clc
            res = input('Press enter to use another program or c to close','s');
            switch res
                case ''
                    
                    Model_program1_1;
                case 'c'
                    return
            end
            
        end 
        
    case '2' %%%%% Test different number of neurons
        
        for lol = 1
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            FileData = load(fullfile(tempdir, 'Fixed_L.mat'));
            Fixed_L = FileData.Fixed_L;
            FileData = load(fullfile(tempdir, 'output_data.mat'));
            Output_Data = FileData.out_train;
            FileData = load(fullfile(tempdir, 'sat_data.mat'));
            Sat_Data = FileData.sat_train;
            FileData = load(fullfile(tempdir, 'lin_data.mat'));
            Lin_Data = FileData.lin_train;
            FileData = load(fullfile(tempdir, 'VD_sat.mat'));
            VD_sat = FileData.VD_sat;
            FileData = load(fullfile(tempdir, 'VD_lin.mat'));
            VD_lin = FileData.VD_lin;
            
            testL = Lin_Data(Lin_Data(:,1) == Fixed_L*1e-6,:);
            testS = Sat_Data(Sat_Data(:,1) == Fixed_L*1e-6,:);
            testO = Output_Data(Output_Data(:,1) == Fixed_L*1e-6,:);

            Out_train = [testO(:,2:5) ; testL(:,2),ones(size(testL,1),1)*VD_lin,testL(:,3:4);testS(:,2),ones(size(testS,1),1)*VD_sat,testS(:,3:4)];
            Out_train(:,end) = log(Out_train(:,end));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalize data  %%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Expression : (target_range*((input_data(:,i)-train_input_min(i))/train_input_range(i)))+target_min; %function minmax
            
            %preprocess data so that minimum is -1 and maximum is 1
            
            [norm_data_train,train_target_min,train_target_range] = norma(Out_train,'define');
            %interval with the number of neurons you want to train the neural network with
            
            max_number_neurons = input('Maximum number of neurons to test?  ]1,+Inf[ _ ') ;
            min_number_neurons = input(['Minimum number of neurons to test?  ]1,',num2str(max_number_neurons),'[ _ ']) ;
            interv_neurons = input('Step size of neurons to test? _ ') ;
            test_neurons = min_number_neurons:interv_neurons:max_number_neurons;
            c = input('Would you like to:\n[1] - test several training percentages \n[2] - perform monte carlo analysis on several number of neurons\n ','s');
            
            switch c
                case '1'
                    
                    %interval with the percentages you want to train the neural network with
                    max_P = input(']0,1[ Maximum training relation to test?   ');
                    min_P = input([']0,',num2str(max_P),'[ Minimum training relation to test?   ']);
                    interv_P = input('Step size to test?  ');
                    P = min_P:interv_P:max_P; %percentage of data training
                case '2'
                    max_P = input('[1,+Inf[ Number of Monte Carlo runs?   ');
                    P = max_P; %percentage of data training
            end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% Train the artificial network  %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Testing severall combinations of number of neurons and differente relation between train and data
        
        search_ann(norm_data_train,test_neurons,P,c) ;
        clc
        
        res = input('Press enter to use another program or c to close','s');
        switch res
            case ''
                Model_program1_1;
            case 'c'
                return
        end
        end   
        
    case '3' %%%%% Create ANN model
        
        for lol = 1
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get User Info %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        ox_t = input('What is the oxide thickness[nm]? ');
        save(fullfile(tempdir, 'ox_t.mat'), 'ox_t', '-mat');
        nneu = input('What is the number of neurons to test? ');
        save(fullfile(tempdir, 'nneu.mat'), 'nneu', '-mat');
        nper = input('What is the training ratio to test? ');
        clc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get stored Data %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileData = load(fullfile(tempdir, 'Fixed_L.mat'));
        Fixed_L = FileData.Fixed_L;
        FileData = load(fullfile(tempdir, 'sat_data.mat'));
        Sat_Data = FileData.sat_train;
        FileData = load(fullfile(tempdir, 'lin_data.mat'));
        Lin_Data = FileData.lin_train;
        FileData = load(fullfile(tempdir, 'output_data.mat'));
        Output_Data = FileData.out_train;
        FileData = load(fullfile(tempdir, 'CV_data.mat'));
        CV_Data = FileData.tr_CV_train;
        FileData = load(fullfile(tempdir, 'VD_sat.mat'));
        VD_sat = FileData.VD_sat;
        FileData = load(fullfile(tempdir, 'VD_lin.mat'));
        VD_lin = FileData.VD_lin;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% Get COV and COX from CV_data file %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% this function returns the Width (um), Cov (pF), Cox (nF/cm^2) from each
        %%% individual transistor size averaging values for COV and Cox
        
        [cap_info] = get_cov(CV_Data,Fixed_L);
        save(fullfile(tempdir, 'cap_info.mat'), 'cap_info', '-mat');        
        average_cox = mean(cap_info(:,3));
        stand_cox = std(cap_info(:,3),0,1);
        
        display(" Cox = " + num2str(round(average_cox,4)) + char(177) + num2str(round(stand_cox,4)) + " nF/cm^2");
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% Get VT from Lin_data file %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% this function returns the Width (um), Cov (pF), Cox (nF/cm^2) from each
        
        
        [vt_info_Lin] = get_vt(Lin_Data,Fixed_L,'lin');
        save(fullfile(tempdir, 'vt_info_Lin.mat'), 'vt_info_Lin', '-mat');        
        average_vt = mean(vt_info_Lin(:,3));
        stand_vt = std(vt_info_Lin(:,3),0,1);
      
        display("Linear VT = " + num2str(round(average_vt,4)) + char(177) + num2str(round(stand_vt,4)) + " V");
        delete (fullfile(tempdir, 'norma_data.mat'))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% Get VT from Sat_data file %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% this function returns the Width (um), Cov (pF), Cox (nF/cm^2) from each
        %%% individual transistor size averaging values for COV and Cox
        
        [vt_info_Sat] = get_vt(Sat_Data,Fixed_L,'sat');
        save(fullfile(tempdir, 'vt_info_Sat.mat'), 'vt_info_Sat', '-mat');        
        average_vt = mean(vt_info_Sat(:,3));
        stand_vt = std(vt_info_Sat(:,3),0,1);
      
        
        display("Saturation VT = " + num2str(round(average_vt,4)) + char(177) + num2str(round(stand_vt,4)) + " V");
        
        delete (fullfile(tempdir, 'norma_data.mat'))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% Get Curve model %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% this function returns the Width (um), Cov (pF), Cox (nF/cm^2) from each
        %%% individual transistor size averaging values for COV and Cox
             
        testL = Lin_Data(Lin_Data(:,1) == Fixed_L*1e-6,:);
        testS = Sat_Data(Sat_Data(:,1) == Fixed_L*1e-6,:);
        testO = Output_Data(Output_Data(:,1) == Fixed_L*1e-6,:);
        
%         Out_train = testO(:,2:5) ;
        
         Out_train = [testO(:,2:5) ; testL(:,2),ones(size(testL,1),1)*VD_lin,testL(:,3:4) ; testS(:,2),ones(size(testS,1),1)*VD_sat,testS(:,3:4)];
         Out_train(:,end) = log(Out_train(:,end));
        mins = min(Out_train);
        range = max(Out_train)-min(Out_train);
        
        save(fullfile(tempdir, 'mins.mat'), 'mins', '-mat');        
        save(fullfile(tempdir, 'range.mat'), 'range', '-mat');        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalize data  %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Expression used: (target_range*((input_data(:,i)-train_input_min(i))/train_input_range(i)))+target_min; %function minmax
        
        %preprocess data so that minimum is -1 and maximum is 1
        
        [norm_data_train,train_target_min,train_target_range] = norma(Out_train,'define');
        
        
        close all % figures. Can be altered to save figures instead?
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% Train the artificial neural network %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [rmse_train,rmse_val,net] = get_ann(norm_data_train,nneu,nper);
        
        display("rmse_train = " + num2str(rmse_train));
        display("rmse_val = " + num2str(rmse_val));
        
      
        save(fullfile(tempdir, 'net.mat'), 'net', '-mat');
        disp('Data measured: ');
        disp(['     Length : ',num2str(Fixed_L),' ',char(181),'m']);
        disp(['     Width : [',num2str(cap_info(:,1)'),'] ',char(181),'m']);
        disp(['     Overlap capacitance : [',num2str(cap_info(:,2)'),'] pF']);
        disp(['     Threshold Voltage : ',num2str(mean(vt_info_Lin(:,3),1)),' V']);
        disp(['     Dieletric Capacitance : ',num2str(mean(cap_info(:,3),1)),' nF/cm^2']);
        disp(['     Dieletric thickness : ',num2str(ox_t),' nm']);
        disp(['     Number of neurons : ',num2str(nneu)]);
        disp('');
        res = input('Press enter to use another program or c to close','s');
        switch res
            case ''
                Model_program1_1;
            case 'c'
                return
        end
        
        end    
        
    case '4' %%%%% Plot Model vs transistors caratheristics curves
        
        for lol = 1
            
        FileData = load(fullfile(tempdir, 'VD_sat.mat'));
        VD_sat = FileData.VD_sat;
        FileData = load(fullfile(tempdir, 'VD_lin.mat'));
        VD_lin = FileData.VD_lin;
        FileData = load(fullfile(tempdir, 'net.mat'));
        net = FileData.net;
        FileData = load(fullfile(tempdir, 'Fixed_L.mat'));
        Fixed_L = FileData.Fixed_L;
        FileData = load(fullfile(tempdir, 'output_data.mat'));
        Output_Data = FileData.out_train;
        FileData = load(fullfile(tempdir, 'sat_data.mat'));
        Sat_Data = FileData.sat_train;
        FileData = load(fullfile(tempdir, 'lin_data.mat'));
        Lin_Data = FileData.lin_train; 
                
        plot_mean_model_std_out(net,Fixed_L,Output_Data);
        plot_mean_model_std_sl(net,Fixed_L,Sat_Data,VD_sat);
        plot_mean_model_std_sl(net,Fixed_L,Lin_Data,VD_lin); 
        clc
        
%         disp('Save generated network?');
%         a = input('[y,n]?','s');
%         
%         if a == 'y'
%             filename = input('Filename?','s');
%           
%             save(fullfile(cd, ['\networks\' filename '.mat']), 'net', '-mat');      
%                   
%         end
        res = input('Press enter to use another program or c to close','s');
        switch res
            case ''
                Model_program1_1;
            case 'c'
                return
        end
        end
        
    case '5' %%%%%% Check model vs transistors transconductance(gm)and output conductance(gds)
        
        for lol = 1
            
        FileData = load(fullfile(tempdir, 'VD_sat.mat'));
        VD_sat = FileData.VD_sat;
        FileData = load(fullfile(tempdir, 'VD_lin.mat'));
        VD_lin = FileData.VD_lin;
        FileData = load(fullfile(tempdir, 'net.mat'));
        net = FileData.net;
        FileData = load(fullfile(tempdir, 'Fixed_L.mat'));
        Fixed_L = FileData.Fixed_L;
        FileData = load(fullfile(tempdir, 'output_data.mat'));
        Output_Data = FileData.out_train;
        FileData = load(fullfile(tempdir, 'sat_data.mat'));
        Sat_Data = FileData.sat_train;
        FileData = load(fullfile(tempdir, 'lin_data.mat'));
        Lin_Data = FileData.lin_train; 
                
        
        get_gm(net,Fixed_L,Sat_Data,VD_sat);
        get_gm(net,Fixed_L,Lin_Data,VD_lin);
        get_gds(net,Fixed_L,Output_Data);

        res = input('Press enter to use another program or c to close','s');
        switch res
            case ''
                Model_program1_1;
            case 'c'
                return
        end 
        end  
       
    case '6' %%%%% Save model and generate verilog A file 
        
        for lol = 1
        FileData = load(fullfile(tempdir, 'ox_t.mat'));
        ox_t = FileData.ox_t;
        FileData = load(fullfile(tempdir, 'Fixed_L.mat'));
        Fixed_L = FileData.Fixed_L;
        FileData = load(fullfile(tempdir, 'net.mat'));
        net = FileData.net;
        FileData = load(fullfile(tempdir, 'nneu.mat'));
        nneu = FileData.nneu;
        FileData = load(fullfile(tempdir, 'vt_info_Lin.mat'));
        vt_info_Lin = FileData.vt_info_Lin;
        FileData = load(fullfile(tempdir, 'cap_info.mat'));
        cap_info = FileData.cap_info;
        FileData = load(fullfile(tempdir, 'mins.mat'));
        mins = FileData.mins;
        FileData = load(fullfile(tempdir, 'range.mat'));
        range = FileData.range;
        
        n_vA_1_1(net,cap_info,vt_info_Lin,Fixed_L,mins,range,ox_t);


        disp('The verilogA file containing the generated model obtained can be found in ');
        disp([cd,'\models']);
        disp('');
        res = input('Press enter to use another program or c to close','s');
        switch res
            case ''
                Model_program1_1;
            case 'c'
                return
        end
        end
        
    case '7'
        
        close all 
        clc
        return
end

if res == 'c'
    return
end



end
