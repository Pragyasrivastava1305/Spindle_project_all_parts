close all 
clear all

load D=0.01_va=0.9_tauc=0.mat;

imagesc(t_reg,x,cL_arr); 
set(gca,'fontsize',14); 
axis xy;  axis equal tight;
colorbar; caxis([0 13]);
xlabel('Time (min)', 'FontSize',16);
ylabel('\bar{c}','FontSize',16,'Interpreter','tex')


hold on; 
plot(t_reg,xN_per_reg,'.w','Markersize',10)
plot(t_reg,xC_per_reg,'.k','Markersize',10)
