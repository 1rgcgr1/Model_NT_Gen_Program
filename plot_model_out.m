function plot_model_out(net,L,train)

close all
format long
clc
% disp('Do you want to plot the mean +- standart deviation and model expected prediction?');
% fig = input('[y,n]? ','s');
y1=zeros(7);
number_std = 1;
p = 1;
s = 1;
disp('Program running')

leng = train(:,4)== train(1,4);

plotVD = train(1:sum(leng)/sum( train(leng,3)==train(1,3)),3);



result_L = (train(:,1)== L*1e-6) ;

Lorg = train(result_L,:); %%% ALL similar Ls together


while any(Lorg)% Separting between Ws
    result_W = (Lorg(:,2)==(Lorg(1,2))) ;
    
    Worg = Lorg(result_W,:);%%% ALL similar Ls,Ws together
    
    f=1;
   
    while any (Worg) % Separating between lines
       
        result_VGS = (Worg(:,4)==(Worg(1,4)));
        VGSorg =  Worg(result_VGS,:);%%% ALL similar Ls,Ws,VGSs together
        k=1;
       
        while any(any(VGSorg))% Separating betwen points
            result_VDS = (VGSorg(:,3)==(VGSorg(1,3)));
            
            VDSorg = VGSorg(result_VDS,:); %%% ALL similar Ls,Ws,VGSs,VDSs together
            M = mean(VDSorg,1);
            S = std(VDSorg(:,5),0,1);
            num_samples = VDSorg(size(VDSorg,2)); %number of samples
            
            if k==1
                y1= [M,S,num_samples];
                
            else
                y1 = [y1;M,S,num_samples];
            end
            
            VGSorg = VGSorg(~result_VDS,:);
            
            k = k+1;
        end
        
        y_call_ann = y1(:,2:4); %uma linha do gráfico que eu vou chamar linha a linha à rede neuronal
        
        if y_call_ann(1,3)>=0
        %%%Ask ANN
        m_line_test = ask_ann(y_call_ann,net);
        %%%    
                
        if f == 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%To correct
            lg = y1(1,4);
            yr = y1(:,5);
            er = y1(:,6);
            mr = m_line_test ;
            %%% asking the neural network vctor
            
        else
            lg = [lg;y1(1,4)]; %legend
            yr = [yr,y1(:,5)]; %original data
            er = [er,y1(:,6)]; %error data = 1 std deviation    
            mr = [mr,m_line_test]; %model data
        end
        end
        Worg = Worg(~result_VGS,:);
        f = f+1;
        
    end
    
    figure
    set(gcf,'color','w');
    VD = plotVD(1:p:end);

    for o = 1:s:size(yr,2)
        
        model_y = mr(1:p:end,o);
%         if model)
        te=plot(VD,model_y,'.');
        te.MarkerSize = 14;
         if o == 1
            fLeg = ["VGS_m = " +  lg(o) + " V"];
        else
            fLeg = [fLeg,"VGS_m = " +  lg(o) + " V"];
        end
        hold on
    end
    
    legend(fLeg,'Location','bestoutside'); %legend for original data
    
    xlabel('V_{DS} (V)')
    ylabel('I_{DS} (A)')
%     ylim([min(min(y)),max(max(y))+number_std*max(max(error))]);
    grid on
    set(gca,'fontname','times','fontsize',12);
    title({"L = " + num2str(round(Lorg(1,1)*1e6)) + " \mum | W = " + num2str(round(Lorg(1,2)*1e6)) + " \mum";"n = " + num2str(num_samples)},'FontName','Times','FontSize',14)
    hold off
    saveas(gcf,[cd,'\Figures\MODEL2_OUT\average_L',num2str(Lorg(1,1)*1e6),'_W',num2str(Lorg(1,2)*1e6),date,'.png']);
    
    Lorg = Lorg(~result_W,:);
    
end


end

