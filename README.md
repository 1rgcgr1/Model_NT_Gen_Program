# Model_NT_Gen_Program
Generate n-type transistor models using neural networks with Matlab 2018R
% 1.1st version
% Author: Rui Guilhoto
% Instructions: 
% Warning: this program does not contain the required code to run on all computers 
% That functionality can be adapted by the user. For now the program can run entirely on command window on my pc
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
% [6] - This step will generate the verilogA file 
%
% [7] - Close the program
%
% For more information
% please contact me by email: r.guilhoto@campus.fct.unl.pt
