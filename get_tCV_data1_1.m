function [ytt] = get_tCV_data1_1
clear all

fi = dir('testes rui\CV_data');
fi = fi(3:end);

disp('%Introduce the number of header lines of one of the files contained in folder');
disp(fi(1).folder);
open(fi(1).name);

%initial line of datavalues in the csv file
ildv = input('? ');

close all
d = csvread(fi(1).name,ildv,1); %reads one file

sd = size(d,1); % number of rows in transistors current file
clc
disp('Do you want to plot the original files data')
fig = input('[y,n]?','s');
clc
disp('Obtaining CV curves');
for i = 1 : size(fi)
    
    [fW,fL,zone] = filename_reader(fi,i,'cov'); % design values in um
    
    data = csvread(fi(i).name,ildv,1);
    
    Cap = data(:,3);
    VG_SD = data(:,1);
    
    switch fig
        case'y'
            figure
            set(gcf,'color','w');
            plot(VG_SD,Cap);
            xlabel('V (V)');
            ylabel('Capacitance (F)');
            ylim([min(min(Cap)),max(max(Cap))]);
            xlim([min(VG_SD),max(VG_SD)]);
            set(gca,'fontname','times','fontsize',12);
            %     legend("C_o_x = "+num2str(Cox)+" nF/cm^2 | C_o_v = " + num2str(Cov)+ " pF",'Location','best');
            %     text((max(VGSD)+min(VGSD))/2,max(max(Cap))+ 0.001*max(max(Cap)),"W = "+ num2str(fW) + " \mum " + "L = "+ num2str(fL) + " \mum" ,'FontName','Times','FontSize', 14,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
            title("L = "+ num2str(fL) + " \mum | W = "+ num2str(fW) + " \mum | zone "+ zone);
    end
    for j = 1:sd
        Twidth(j,1) = fW *1e-6;
        Tlength(j,1) = fL*1e-6;
    end
    % Store information
    if i == 1
        store_cap = [Tlength,Twidth,VG_SD,Cap];
    else
        store_cap = [store_cap;Tlength,Twidth,VG_SD,Cap];
    end
end
clc

disp('Do you want to plot the mean +- standart deviation');
fig = input('[y,n]? ','s');
close all
number_std = 1;
p = 10;
s = 3;
disp('Program running')
i=1;
while any(store_cap)%separating between Ls
    
    result_L = (store_cap(:,1)==(store_cap(1,1)));
    
    Lorg = store_cap(result_L,:); %%% ALL similar Ls together
    
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
                errorbar(x,y,error);
                
                xlabel('VG_S_D (V)')
                ylabel('Capacitance (C)')
                ylim([min(min(y))-number_std*max(max(error)),max(max(y))+number_std*max(max(error))]);
                set(gca,'fontname','times','fontsize',12);
                title({"L = " + num2str(round(Lorg(1,1)*1e6)) + " \mum | W = " + num2str(round(Lorg(1,2)*1e6)) + " \mum";"n = " + num2str(num_samples)},'FontName','Times','FontSize',14)
                hold off
                saveas(gcf,[cd,'\Figures\CV\average_L',num2str(Lorg(1,1)*1e6),'_W',num2str(Lorg(1,2)*1e6),date,'.png']);
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
    store_cap = store_cap(~result_L,:);
    i = i+1;
end
end