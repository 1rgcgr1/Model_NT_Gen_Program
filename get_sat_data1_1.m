function [ytt,VD] = get_sat_data1_1

clear all

fi = dir('testes rui/Saturation_data');
fi = fi(3:size(fi));



disp('%Introduce the number of header lines of the files contained in folder'); 
disp(fi(1).folder);
open(fi(1).name);
%initial line of datavalues in the csv file 
ildv = input('? ');

close all

d = csvread(fi(1).name,ildv,1); %reads one file 
 
ppl = 1;
 

% VGS_step_size   = round(d(ppl+1,2)-d(ppl,2),3); % step size of VGS 
% VDS_step_size   = round(d(2,1)-d(1,1),3); % step size of VDS 

% VDS_values = [d(1,2):VDS_step_size:d(end,2)]; %define minimum and maximum VD values
% VGS_values = [d(1,1):VGS_step_size:d(end,1)]; %define minimum and maximum VG values
total_data = zeros(1,1);
VD = d(1,2);
% number_std = 1;
sd = size(d,1); % number of rows in transistors current file
clc
disp('Do you want to plot the original files data');
fig = input('[y,n]?','s');
clc
disp('Obtaining Saturation curves');
for i = 1:size(fi)
    
    [fW,fL,zone] = filename_reader(fi,i,'sat'); % gets the transistors sizes in filename
     
    d = csvread(fi(i).name,ildv,1);
%     VGS_step_size   = round(d(ppl+1,2)-d(ppl,2),3);
   
    VG = d(:,1);
    ID = d(:,8);

    
    switch fig
        case'y'
            figure
            set(gcf,'color','w');
            plot(VG,ID)
            xlabel('VDS (V)')
            ylabel('I_{DS} (A)')
            ax = gca;
            set(ax,'fontname','times','fontsize',12);
            ax.YScale = 'log';
            grid on
            xticks([-10:2:10]);
            ylim([5e-14,1e-4]);
            title({"L = " + num2str(fL) + " \mum | W = " + num2str(fW) + " \mum";"zone " + zone},'FontName','Times','FontSize',14)
            
    end
    for j = 1:sd
        Twidth(j,1) = fW *1e-6;
        Tlength(j,1) = fL *1e-6;
    end
    %Storing information in datastore matrix
    data = [ Tlength ,Twidth, VG, ID ]; %current is in uA
    
    if i == 1
        total_data = data ;
    else
        total_data = [total_data;data];
    end

end
clc

disp('Do you want to plot the mean +- standart deviation');
fig = input('[y,n]? ','s');
close all
number_std = 1;
p = 10;
s = 3;

%  sat_train = [ Tlength ,Twidth, VG, AbsID ]; %current is in uA
disp('Program running')
i=1;
while any(total_data)%separating between Ls
    
    result_L = (total_data(:,1)==(total_data(1,1)));
    
    Lorg = total_data(result_L,:); %%% ALL similar Ls together
    
    j=1;
    while any(Lorg)% Separting between Ws
        
        result_W = (Lorg(:,2)==(Lorg(1,2))) ;
        
        Worg = Lorg(result_W,:);%%% ALL similar Ls,Ws together
        
        f=1;
        
        
        while any (Worg) % Separating between lines
            
            result_VGS = (Worg(:,3)==(Worg(1,3)));
            
            VGSorg =  Worg(result_VGS,:);%%% ALL similar Ls,Ws,VGSs together
            
            M = mean(VGSorg,1);
            S = std(VGSorg(:,4),0,1);
            num_samples = size(VGSorg,1); %number of samples
            
            if f==1
                y1 = [M,S,num_samples];
                % the mean, standart deviation and nunber of samples in Id for each value of L,W,VG
                
            else
                y1 = [y1;M,S,num_samples];
            end
            
            
            
            Worg = Worg(~result_VGS,:);
            f = f+1;
            
        end
        
        switch fig
            case 'y'
                
                figure
                set(gcf,'color','w');
                error = number_std*y1(1:p:end,5);
                y = y1(1:p:end,4);
                x = y1(1:p:end,3);
                af = errorbar(x,y,error);
                af.LData = af.YData - max(1e-13,af.YData-af.LData); % eliminates negative values by comparing floating point accuracy with errorbar values
                xlabel('VGS (V)')
                ylabel('I_{DS} (A)')
                ax = gca;
                set(ax,'fontname','times','fontsize',12);
                ax.YScale = 'log';
                grid on
                xticks([-10:2:10]);
                ylim([5e-14,1e-4]);
                title({"L = " + num2str(Lorg(1,1)*1e6) + " \mum | W = " + num2str(Lorg(1,2)*1e6) + " \mum";"n = " + num2str(num_samples)},'FontName','Times','FontSize',14)
                hold off
                saveas(gcf,[cd,'\Figures\SAT\average_L',num2str(Lorg(1,1)*1e6),'_W',num2str(Lorg(1,2)*1e6),date,'.png']);
               
        end
        if j==1
            yt = y1;
            % the mean, standart deviation and nunber of samples in Id for each value of L,W,VG
        else
            yt = [yt;y1];
        end
        
        Lorg = Lorg(~result_W,:);
        j = j+1;
    end
    
    if i==1
            ytt = yt;
            % the mean, standart deviation and nunber of samples in Id for each value of L,W,VG
        else
            ytt = [ytt;yt];
        end
    
    total_data = total_data(~result_L,:);
    i = i+1;
end
end