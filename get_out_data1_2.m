function [ytt] = get_out_data1_2
format long
clear all

fi = dir('testes rui/Output_data');
fi = fi(3:size(fi));



disp('%Introduce the number of header lines of the files contained in folder');
disp(fi(1).folder);
open(fi(1).name);
%initial line of datavalues in the csv file
ildv = input('? ');

close all

d = csvread(fi(1).name,ildv,1); %reads one file

ppl = 1;

while d(ppl,2) == d(ppl+1,2)
    ppl = ppl+1;
end

VGS_step_size   = round(d(ppl+1,2)-d(ppl,2),3); % step size of VGS
VDS_step_size   = round(d(2,1)-d(1,1),3); % step size of VDS

VDS_values = [d(1,2):VDS_step_size:d(end,2)]; %define minimum and maximum VD values
VGS_values = [d(1,2):VGS_step_size:d(end,2)]; %define minimum and maximum VG values
total_data = zeros(1,1);

sd = size(d,1); % number of rows in transistors current file
clc
disp('Do you want to plot the original files data')
fig = input('[y,n]?','s');
clc
disp('Obtaining Output curves');
for i = 1:size(fi)
    
    [fW,fL,zone] = filename_reader(fi,i,'out'); % gets the transistors sizes in filename
    
    d = csvread(fi(i).name,ildv,1);
    
    VD = d(:,1);
    VG = d(:,2);
    AbsID = d(:,7);
    plotVD = VD(1:ppl,1);
    ncurves = size(AbsID,1)/ppl;
    
    for o = 0:ncurves-1
        if o == 0
            IDcurves = AbsID(1:ppl,1);
            legs = ["VGS = "+VGS_values(1)+" V"];
        else
            IDcurves = [IDcurves,AbsID((o*ppl)+1:(1+o)*ppl,1)];
            legs = [legs,"VGS = "+VGS_values(1+o)+" V"];
        end
    end
    
    switch fig
        case'y'
            figure
            set(gcf,'color','w');
            plot(plotVD,IDcurves)
            legend(legs,'Location','best');
            xlabel('VDS (V)')
            ylabel('I_{DS} (A)')
            ylim([min(min(IDcurves)),max(max(IDcurves))]);
            set(gca,'fontname','times','fontsize',12);
            text(max(plotVD)/2,max(max(IDcurves))+ 0.07*max(max(IDcurves)),"Training IV characteristics",'FontName','Times','FontSize', 14,'FontWeight','Bold','HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
            text(max(plotVD)/2,max(max(IDcurves))+ 0.005*max(max(IDcurves)),"L = " + num2str(fL) + " \mum | W = " + num2str(fW) + " \mum | zone " + zone,'FontName','Times','FontSize', 14,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
            set(gca,'InnerPosition',[0.13 0.11 0.775 0.77])
            
    end
    for j = 1:sd
        Twidth(j,1) = fW *1e-6;
        Tlength(j,1) = fL*1e-6;
    end
    %Storing information in datastore matrix
    data = [ Tlength ,Twidth, VD, VG, AbsID ]; %current is in uA
    
    if i == 1
        total_data = data ;
    else
        total_data = [total_data;data];
    end
    %  lin_train = [ Tlength ,Twidth, VD, VG, AbsID ]; %current is in uA
end 

clc
data = total_data;
disp('Do you want to plot the mean +- standard deviation');
fig = input('[y,n]? ','s');
disp('Program running')
close all
number_std = 1;
p = 10;
s = 3;
i=1;

plotVD = total_data(1:(size(total_data(:,3),1)/sum( total_data(:,3)==total_data(1,3) )),3);
while any(total_data)%separating between Ls
    
    result_L = total_data(:,1)==(total_data(1,1));
    
    Lorg = total_data(result_L,:); %%% ALL similar Ls together
    
    j=1;
    while any(Lorg)% Separting between Ws
        result_W = (Lorg(:,2)==(Lorg(1,2))) ;
        
        Worg = Lorg(result_W,:);%%% ALL similar Ls,Ws together
        
        f=1;
        
        while any (Worg) % Separating between lines
            
            result_VGS = (Worg(:,4)==(Worg(1,4)));
            VGSorg =  Worg(result_VGS,:);%%% ALL similar Ls,Ws,VGSs together
            k=1;
            
            while any(VGSorg(:,3))% Separating betwwen points
                result_VDS = (VGSorg(:,3)==(VGSorg(1,3)));
                
                VDSorg = VGSorg(result_VDS,:); %%% ALL similar Ls,Ws,VGSs,VDSs together
                M = mean(VDSorg,1);
                S = std(VDSorg(:,5),0,1);
                num_samples = size(VDSorg,1); %number of samples
                
                if k==1
                    y1= [M,S,num_samples];
                    
                else
                    y1 = [y1;M,S,num_samples];
                end
                
                VGSorg = VGSorg(~result_VDS,:);
                
                k = k+1;
            end
            
            if f == 1
                lg = y1(1,4);
                yr = y1(:,5);
                er = y1(:,6);
                ya = y1;
            else
                lg = [lg;y1(1,4)];
                yr = [yr,y1(:,5)];
                er = [er,y1(:,6)];
                ya = [ya;y1];
            end
            
            Worg = Worg(~result_VGS,:);
            f = f+1;
            
        end
        switch fig
            case 'y'
                figure
                set(gcf,'color','w');
                for o = 1:s:size(yr,2)
                    error = number_std*er(1:p:end,o);
                    y = yr(1:p:end,o);
                    VD = plotVD(1:p:end);
                    errorbar(VD,y,error);
                    hold on
                    if o == 1
                        fLeg = ["VGS_d = " +  lg(o) + " V"];
                    else
                        fLeg = [fLeg,"VGS_d = " +  lg(o) + " V"];
                    end
                    
                end
                
                legend(fLeg,'Location','best'); %legend for original data
                
                xlabel('VDS (V)')
                ylabel('I_{DS} (A)')
                ylim([min(min(y)),max(max(y))+number_std*max(max(error))]);
                set(gca,'fontname','times','fontsize',12);
                title({"L = " + num2str(round(Lorg(1,1)*1e6)) + " \mum | W = " + num2str(round(Lorg(1,2)*1e6)) + " \mum";"n = " + num2str(num_samples)},'FontName','Times','FontSize',14)
                hold off
                saveas(gcf,[cd,'\Figures\OUT\average_L',num2str(round(Lorg(1,1)*1e6)),'_W',num2str(round(Lorg(1,2)*1e6)),date,'.png']);
        end
        if j == 1
            yt = ya;
        else
            yt = [yt;ya];
        end
        Lorg = Lorg(~result_W,:);
        j=j+1;
    end
    
    if i == 1
        ytt = yt;
    else
        ytt = [ytt;yt];
    end
    total_data = total_data(~result_L,:);
    i=i+1;
end
end
