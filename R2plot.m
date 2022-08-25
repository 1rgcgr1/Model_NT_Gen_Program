function  R2plot(sampleyTrain,modelyTrain,sampleyVal,modelyVal,R2_t,R2_v,h)

figure
set(gcf,'color','w');
plot(sampleyTrain,modelyTrain,'x')
hold on
plot(-1:1,-1:1)
hold on
plot(sampleyVal,modelyVal,'o')
ylim([-1 1]);
xlim([-1 1]);
ylabel('Normalized predicted data')
xlabel('Normalized training data')
text(0,1.20,'Correlation between training data and model data','FontName','Times','FontSize', 14,'FontWeight','Bold','HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
text(0,1.00,"R^2_t = " + num2str(R2_t(h)) +  " | R^2_v =  " + num2str(R2_v(h)) ,'FontName','Times','FontSize', 13,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
set(gca,'InnerPosition',[0.13 0.11 0.775 0.76])
set(gca,'FontName','Times','FontSize',11)
hold off

end