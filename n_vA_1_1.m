function n_vA_1_1(net,cap_info,vt_info,FixedL,mins,range,ox_t)

clc
disp('Saving model in txt file');
disp('Set a model name. Cannot start with a number!');
filename_model = input('','s');

format shorteng

number_of_layers = size(net.b,1);

wts = [net.IW,net.LW];  %  all layer weights
biases = net.b; %all layer bias

hidden_weights = wts(1:number_of_layers-1,1:number_of_layers-1);  %size NNHL x 1 hidden layer weights
hidden_bias = biases(1:number_of_layers-1);        %size NNHL x 1 hidden layer bias

out_weights = wts(number_of_layers,number_of_layers);    %size 1 x NNHL output layer weights
out_bias = biases(number_of_layers);           % size 1 x 1 output layer bias


w = cap_info(:,1); 
l = cap_info(:,2);
cov = cap_info(:,2);
cox = mean(cap_info(:,3));

all_vt = vt_info(:,3);
cd('models');
fid = fopen( [filename_model ,'_',date, '.txt'], 'w');
i=1;
g=1;
j=1;
p=1;

% declaring initial includes and ports

fprintf(fid,'`include "constants.h" ');
fprintf(fid,'\n`include "disciplines.h" ');
fprintf(fid,['\n\nmodule ',filename_model,'(D,G,S); ']);
fprintf(fid,' \ninout G,D,S ;');
fprintf(fid,' \nelectrical D,G,S;');

% declare parameterswidth,legth,dieletric thickness,COX,vth 

fprintf(fid,'\nparameter integer NI = 3;');
for i=1:number_of_layers-1
fprintf(fid,['\nparameter integer NNHL',num2str(i),' = %d; '],size(net.b{i},1));
end
fprintf(fid,'\nparameter integer NO = 1 ;');
fprintf(fid,'\nparameter real Width = %.9f from (0:inf);',FixedL*1e-6);
fprintf(fid,'\nparameter real Length = %.9f from (0:inf);',FixedL*1e-6);
fprintf(fid,'\nparameter real tickness = %.12f from (0:inf);',ox_t*1e-9); 
fprintf(fid,'\nparameter real epsila = %.15f from (0:inf);',8.851*1e-12); % fixed variable
fprintf(fid,'\nparameter real COX = %.8f from (0:inf);',cox*1e-5);
fprintf(fid,'\nparameter real vth = %.3f from (-5:inf);',mean(all_vt));

% declare hidden layers weights and biases

if number_of_layers == 2
    
        g=1;
        fprintf(fid,['\nreal hlayer1_w[0:(NI*NNHL1)-1];']);
        fprintf(fid,'\nreal hlayer1_b[0:NNHL1-1];');
        fprintf(fid,'\nreal hlayer1_y[0:NNHL1-1];');
        fprintf(fid,'\nreal hlayer1_v[0:NNHL1-1];');
        
        fprintf(fid,['\nreal olayer_w[0:(NO*NNHL',num2str(g),')-1];']);
        fprintf(fid,'\nreal olayer_b[0:NO-1];');
        fprintf(fid,'\nreal olayer_y[0:NO-1];');
        fprintf(fid,'\nreal olayer_v[0:NO-1];');
        
else     
        
        fprintf(fid,['\nreal hlayer1_w[0:(NI*NNHL1)-1];']);
        fprintf(fid,'\nreal hlayer1_b[0:NNHL1-1];');
        fprintf(fid,'\nreal hlayer1_y[0:NNHL1-1];');
        fprintf(fid,'\nreal hlayer1_v[0:NNHL1-1];');
        
        for g = 2:number_of_layers-1
            fprintf(fid,['\nreal hlayer',num2str(g),'_w[0:(NNHL',num2str(g-1),'*NNHL',num2str(g),')-1];']);
            fprintf(fid,['\nreal hlayer',num2str(g),'_b[0:NNHL',num2str(g),'-1];']);
            fprintf(fid,['\nreal hlayer',num2str(g),'_y[0:NNHL',num2str(g),'-1];']);
            fprintf(fid,['\nreal hlayer',num2str(g),'_v[0:NNHL',num2str(g),'-1];']);
        end
        
        fprintf(fid,['\nreal olayer_w[0:(NO*NNHL',num2str(g),')-1];']);
        fprintf(fid,'\nreal olayer_b[0:NO-1];');
        fprintf(fid,'\nreal olayer_y[0:NO-1];');
        fprintf(fid,'\nreal olayer_v[0:NO-1];');
end
%declare variables

fprintf(fid,'\nreal width;');
fprintf(fid,'\nreal cgs,cgd, qgs, qgd,cch_d, cch_semi, cch, cov_d, cov_semi, cov, vgs,vds;');
fprintf(fid,'\ninteger i,j,ii,jj,iii,jjj,k;');

%declare pre and post processing parameters

fprintf(fid,'\nreal train_input_range[0:(NI-1)] ;');
fprintf(fid,'\nreal train_input_min[0:(NI-1)] ;');
fprintf(fid,'\nreal train_output_range,id;');
fprintf(fid,'\nreal train_output_min;');
 

%analog begin

fprintf(fid,'\n\nanalog begin');
fprintf(fid,'\n\n\t@(initial_step or initial_step("static"))');
fprintf(fid,'\n\n\nbegin\n\n');

%overlap capacitances

fprintf(fid,'\n\n\tif(Width == %f)',w(1)*1e-6);
fprintf(fid,'\nbegin');
fprintf(fid,'\n\tcov = %.15f ; ',cov(1)*1e-12);
fprintf(fid,'\nend');

for i = 2 : size(cap_info,1)
    fprintf(fid,'\n\n\telse if(Width == %f)',w(i)*1e-6);
    fprintf(fid,'\nbegin');
    fprintf(fid,'\n\tcov = %.15f ; ',cov(i)*1e-12);
    fprintf(fid,'\nend');
end

fprintf(fid,'\n\n\tend ');

fprintf(fid,'\n\n cch = COX*Width*Length;');

%modelnetwork hidden weights

fprintf(fid, '\n\n //hidden weights ');


for j=1:number_of_layers-1
    fprintf(fid, ['\n //layer ',num2str(j)]);
    number_of_neuron_per_h_layer = size(hidden_weights{j,j},1);
    number_of_inputs_per_h_layer = size(hidden_weights{j,j},2);
    value_hw = hidden_weights{j,j};
    k=0;
    
    for i=1:number_of_neuron_per_h_layer
        fprintf(fid, '\n');
        
        for p = 1: number_of_inputs_per_h_layer
            fprintf(fid, ['hlayer',num2str(j),'_w[ %d] =  %f;'], k , value_hw(i,p));
            fprintf(fid, '\t');
            k=k+1;
        end
        
    end
    fprintf(fid, '\n');
end
%modelnetwork hidden bias
fprintf(fid, '\n\n\n //hidden bias ');

for j=1:number_of_layers-1
    fprintf(fid, ['\n\n //layer ',num2str(j)]);
    number_of_neuron_per_h_layer = size(hidden_bias{j,1},1);
    value_b = hidden_bias{j,1};
    k=0;
    for i=1:number_of_neuron_per_h_layer
        fprintf(fid, ['\nhlayer',num2str(j),'_b[ %d] =  %f;'],k, value_b(i));
        k=k+1;
    end
end
%modelnetwork output weights

fprintf(fid, '\n\n //output layer weights ');
fprintf(fid, '\n');

value_weigths_output = out_weights{1,1};
for i = 1:size(value_weigths_output,2)
   fprintf(fid, '\nolayer_w[ %d] =  %f; ',i-1, value_weigths_output(i));
end

%modelnetwork output bias

fprintf(fid, '\n\n\n //output layer bias ');
fprintf(fid, '\n\n\n');


fprintf(fid, 'olayer_b[ 0] =  %f ; \n', out_bias{1,1});

fprintf(fid, '\n\n\n');

%Declaring ranges for (de)normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(range,2)-1
fprintf(fid, '\ntrain_input_range[%.0f] = %.9f;',i-1, range(i));
end

fprintf(fid, '\ntrain_output_range = %f;',range(length(range)));

%Declaring minimums for (de)normalization

for i=1:size(mins,2)-1
fprintf(fid, '\ntrain_input_min[%.0f] = %.9f;',i-1, mins(i));
end

fprintf(fid, '\ntrain_output_min = %.15f;',mins(length(mins)));

%Declaring hiddenlayer variables values 

for i= 1:number_of_layers-1
fprintf(fid, ['\n\nfor (i = 0; i < (NNHL',num2str(i),'); i = i+1)']);
fprintf(fid, '\nbegin');
fprintf(fid, ['\n\thlayer',num2str(i),'_v[i] = 0;']);
fprintf(fid, '\nend');
end

fprintf(fid, '\n\nfor (i = 0; i < (NO); i = i+1)');
fprintf(fid, '\nbegin');
fprintf(fid, '\n\tolayer_v[i] = 0;');
fprintf(fid, '\nend');

%Normalizing user defined data width, vds, vgs

fprintf(fid,'\n\nwidth= (2*(Width - train_input_min[0])/train_input_range[0]) -1;');
fprintf(fid,'\nvds    = (2*(V(D,S) - train_input_min[1])/train_input_range[1]) -1;');
fprintf(fid,'\nvgs    = (2*(V(G,S) - train_input_min[2])/train_input_range[2]) -1;');

%Meyer's approximation 

fprintf(fid,'\nif (V(G,S) <  vth) ');
fprintf(fid,'\nbegin');
fprintf(fid,'\n\tcgs = cov;');
fprintf(fid,'\n\tcgd = cov;');
fprintf(fid,'\nend');
fprintf(fid,'\nelse if(V(D,S) >= V(G,S) - vth) ');
fprintf(fid,'\nbegin');
fprintf(fid,'\n\tcgs  = (2.0/3.0)*cch  + cov;');
fprintf(fid,'\n\tcgd  = cov;');
fprintf(fid,'\nend');
fprintf(fid,'\nelse');
fprintf(fid,'\nbegin');
fprintf(fid,'\n\tcgs  = 0.5*cch + cov;');
fprintf(fid,'\n\tcgd  = 0.5*cch + cov;');
fprintf(fid,'\n\nend');
fprintf(fid,'\n\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////');
% inputing the user defined data into the ANN
if number_of_layers == 2
    
        fprintf(fid,'\n\nfor (i = 0; i< (NNHL1); i = i+1)');
        fprintf(fid,'\nbegin');
        fprintf(fid,'\n\thlayer1_v[i] = 0;');
        fprintf(fid,'\n\thlayer1_v[i]  =  hlayer1_v[i] + width*hlayer1_w[NI*(i)] + vds*hlayer1_w[NI*(i)+1] +  vgs*hlayer1_w[NI*(i)+2];');
        fprintf(fid,'\n\thlayer1_v[i]  =  hlayer1_v[i] + hlayer1_b[i];');
        fprintf(fid,'\n\thlayer1_y[i]  =  tanh(hlayer1_v[i]);');
        fprintf(fid,'\nend');
        fprintf(fid,'\n\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////');
        fprintf(fid,'\n\nfor (iii = 0; iii< (NO ); iii = iii+1)');
        fprintf(fid,'\nbegin');
        fprintf(fid,'\n\tolayer_v[iii] = 0;');
        fprintf(fid,['\n\tfor (jjj = 0; jjj< (NNHL',num2str(j),'); jjj = jjj+1)']);
        fprintf(fid,'\n\tbegin');
        fprintf(fid,['\n\t\tolayer_v[iii]  =  olayer_v[iii] +hlayer',num2str(j),'_y[jjj] *olayer_w[jjj] ;']);
        fprintf(fid,'\n\tend');
        fprintf(fid,'\n\tolayer_v[iii]  =  olayer_v[iii]  +  olayer_b[iii] ;');
        fprintf(fid,'\n\tolayer_y[iii]  =  (olayer_v[iii]);');
        fprintf(fid,'\n\tolayer_y[iii]  =  (((olayer_y[iii] +1)/2)*train_output_range) + train_output_min;');
        fprintf(fid,'\n\tid  =  exp(olayer_y[iii]);');
        fprintf(fid,'\nend');
        fprintf(fid,'\n\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////');

else
        % inputing the user defined data into the ANN
        fprintf(fid,'\n\nfor (i = 0; i< (NNHL1); i = i+1)');
        fprintf(fid,'\nbegin');
        fprintf(fid,'\n\thlayer1_v[i] = 0;');
        fprintf(fid,'\n\thlayer1_v[i]  =  hlayer1_v[i] + width*hlayer1_w[NI*(i)] + vds*hlayer1_w[NI*(i)+1] +  vgs*hlayer1_w[NI*(i)+2];');
        fprintf(fid,'\n\thlayer1_v[i]  =  hlayer1_v[i] + hlayer1_b[i];');
        fprintf(fid,'\n\thlayer1_y[i]  =  tanh(hlayer1_v[i]);');
        fprintf(fid,'\nend');
        fprintf(fid,'\n\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////');

        % inputing the user defined data into the ANN 2nd+ layer 
        for j=2:number_of_layers-1
            fprintf(fid,'\nk = 0;');
            fprintf(fid,['\n\nfor (ii = 0; ii< (NNHL',num2str(j),' ); ii = ii+1)']);
            fprintf(fid,'\nbegin');
            fprintf(fid,['\n\thlayer',num2str(j),'_v[ii] = 0;']);
            fprintf(fid,['\n\tfor (jj = 0; jj< (NNHL',num2str(j-1),'); jj = jj+1)']);
            fprintf(fid,'\n\tbegin');
            fprintf(fid,['\n\t\thlayer',num2str(j),'_v[ii]  =  hlayer',num2str(j),'_v[ii] +hlayer',num2str(j-1),'_y[jj] *hlayer',num2str(j),'_w[(k)] ;']);
            fprintf(fid,'\n\t\tk = k+1;');
            fprintf(fid,'\n\tend');
            fprintf(fid,['\n\thlayer',num2str(j),'_v[ii]  =  hlayer',num2str(j),'_v[ii]  +  hlayer',num2str(j),'_b[ii] ;']);
            fprintf(fid,['\n\thlayer',num2str(j),'_y[ii]  =  tanh(hlayer',num2str(j),'_v[ii]);']);
            fprintf(fid,'\nend');
            
            fprintf(fid,'\n\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////');

        end
        % inputing the user defined data into the ANN last layer and denormalization
        fprintf(fid,'\n\nfor (iii = 0; iii< (NO ); iii = iii+1)');
        fprintf(fid,'\nbegin');
        fprintf(fid,'\n\tolayer_v[iii] = 0;');
        fprintf(fid,['\n\tfor (jjj = 0; jjj< (NNHL',num2str(j),'); jjj = jjj+1)']);
        fprintf(fid,'\n\tbegin');
        fprintf(fid,['\n\t\tolayer_v[iii]  =  olayer_v[iii] +hlayer',num2str(j),'_y[jjj] *olayer_w[jjj] ;']);
        fprintf(fid,'\n\tend');
        fprintf(fid,'\n\tolayer_v[iii]  =  olayer_v[iii]  +  olayer_b[iii] ;');
        fprintf(fid,'\n\tolayer_y[iii]  =  (olayer_v[iii]);');
        fprintf(fid,'\n\tolayer_y[iii]  =  (((olayer_y[iii] +1)/2)*train_output_range) + train_output_min;');
        fprintf(fid,'\n\tid  =  exp(olayer_y[iii]);');
        fprintf(fid,'\nend');
        fprintf(fid,'\n\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////');

        
end
% Port attribution

fprintf(fid,'\n\nI(D,S)     <+  	(id);');
fprintf(fid,'\nqgs = cgs*V(G,S);');
fprintf(fid,'\nqgd = cgd*V(G,D);');
fprintf(fid,'\nI(G,S)   <+  	ddt(qgs);');
fprintf(fid,'\nI(G,D)   <+  	ddt(qgd);');
fprintf(fid,'\nend');
fprintf(fid,'\nendmodule');

fclose(fid);

cd ..
end


