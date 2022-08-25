function plotterror(x,y,fitted,rmse_train,rmse_val,c)

switch c
    case '1' 
figure
surf(x,y,fitted);
shading interp
title('Relation between error of training and validation');
ylim('auto');
xlim([min(x),max(x)]);
xticks('auto');
xlabel('Number of neurons');
ylabel('Training data relation');
zlabel('RMSE_t / RMSE_v');
set(gca,'fontname','times','fontsize',12);
saveas(gcf,[cd,'\Figures\MODEL_MULTI_TEST\fitted.png']);
figure
surf(x,y,rmse_train);
shading interp
ylim('auto');
title('Root mean square error for training data')
xlim([min(x),max(x)]);
xticks('auto');
xlabel('Number of neurons');
ylabel('Training data relation');
zlabel('Root mean square error');
set(gca,'fontname','times','fontsize',12);
saveas(gcf,[cd,'\Figures\MODEL_MULTI_TEST\train.png']);

figure
surf(x,y,rmse_val);
shading interp
title('Root mean square error for validation data')
xlim([min(x),max(x)]);
xticks('auto');
xlabel('Number of neurons');
ylim('auto');
ylabel('Training data relation');
zlabel('Root mean square error');
set(gca,'fontname','times','fontsize',12);
saveas(gcf,[cd,'\Figures\MODEL_MULTI_TEST\val.png']);
    case '2'
        
figure
plot(x,fitted);
title('Mean Relation between error of training and validation');
ylim('auto');
xlim([min(x),max(x)]);
xticks('auto');
xlabel('Number of neurons');
ylabel('RMSE_t / RMSE_v');
set(gca,'fontname','times','fontsize',12);
saveas(gcf,[cd,'\Figures\MODEL_MULTI_TEST\mfitted.png']);
figure
plot(x,rmse_train);
ylim('auto');
title('Mean Root mean square error for training data')
xlim([min(x),max(x)]);
xticks('auto');
xlabel('Number of neurons');
ylabel('Root mean square error');
set(gca,'fontname','times','fontsize',12);
saveas(gcf,[cd,'\Figures\MODEL_MULTI_TEST\mtrain.png']);
figure
plot(x,rmse_val);
title('Mean root mean square error for validation data')
xlim([min(x),max(x)]);
xticks('auto');
xlabel('Number of neurons');
ylim('auto');
ylabel('Root mean square error');
set(gca,'fontname','times','fontsize',12);
saveas(gcf,[cd,'\Figures\MODEL_MULTI_TEST\mval.png']);
case '3'
        
figure
plot(x,fitted);
title('Mean Relation between error of training and validation');
ylim('auto');
xlim([min(x),max(x)]);
xticks('auto');
xlabel('Number of neurons');
ylabel('RMSE_t / RMSE_v');
set(gca,'fontname','times','fontsize',12);
saveas(gcf,[cd,'\Figures\MODEL_MULTI_TEST\mfitted2.png']);
figure
plot(x,rmse_train);
ylim('auto');
title('Mean Root mean square error for training data')
xlim([min(x),max(x)]);
xticks('auto');
xlabel('Number of neurons');
ylabel('Root mean square error');
set(gca,'fontname','times','fontsize',12);
saveas(gcf,[cd,'\Figures\MODEL_MULTI_TEST\mtrain2.png']);
figure
plot(x,rmse_val);
title('Mean root mean square error for validation data')
xlim([min(x),max(x)]);
xticks('auto');
xlabel('Number of neurons');
ylim('auto');
ylabel('Root mean square error');
set(gca,'fontname','times','fontsize',12);
saveas(gcf,[cd,'\Figures\MODEL_MULTI_TEST\mval2.png']);

end
end