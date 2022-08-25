function get_gm(net,L,train,VD)


format long
clc

number_std = 1;
p = 1;
s = 4;
%  sat_train = [ Tlength ,Twidth, VG, AbsID ]; %current is in uA
disp('Program running')



result_L = (train(:,1)== L*1e-6 );

Lorg = train(result_L,:); %%% ALL similar Ls together

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
  
    error = number_std*y1(1:p:end,5);
    y = diff(y1(1:p:end,4));
    x = y1(1:p:end,3);
    
    VDm = ones(size(x))*VD;
    Wm = ones(size(x))*Lorg(1,2);
    
    y_call_ann =[Wm,VDm,x];
    
    %%%Ask ANN
    m_line_test = diff(ask_ann(y_call_ann,net));
    %%%    
        
    figure
    set(gcf,'color','w');
    
    af = plot(x(1:size(x,1)-1,:),y);
    xlabel('V_{GS} (V)')
    ylabel('g_m (S)')
    ax = gca;
    set(ax,'fontname','times','fontsize',12);
    
    title({"L = " + num2str(Lorg(1,1)*1e6) + " \mum | W = " + num2str(Lorg(1,2)*1e6) + " \mum"},'FontName','Times','FontSize',14)
    leg = "VDS_d = " +  num2str(VD) + " V";
    
    hold on
    
    te=plot(x(1:size(x,1)-1,:),m_line_test,'.');
    te.MarkerSize = 14;
    leg = [leg, "VDS_m = " +  num2str(VD) + " V"];
    
    
    
    legend (leg,'Location','best');
    ax.YScale = 'linear';
    hold off
    saveas(gcf,[cd,'\Figures\MODEL_GM\average_L',num2str(Lorg(1,1)*1e6),'_W',num2str(Lorg(1,2)*1e6),date,'.png']);

    
    Lorg = Lorg(~result_W,:);
    j = j+1;
end

end
