close all
Dvec=[0 0.01 0.1 1]; 

I0 =  find(param_small(:,1)==0);
D0_dat = param_small(I0,:); 

I1 =  find(param_small(:,1)==0.01);
D001_dat = param_small(I1,:); 

I2 =  find(param_small(:,1)==0.1);
D01_dat = param_small(I2,:); 

I3 =  find(param_small(:,1)==1);
D1_dat = param_small(I3,:); 


figure(1)
subplot(2,2,1)
for i=1:size(I0,1)
    plot(1:3,D0_dat(i,3:5),'o','linewidth',2,'MarkerSize',6); drawnow
    set(gca,'yscale','log')
    set(gca,'xticklabels',[])
    axis square
    xlim([0 4])
    hold on 
end

subplot(2,2,2)
for i=1:size(I1,1)
    plot(1:3,D001_dat(i,3:5),'o','linewidth',2,'MarkerSize',6); drawnow
    set(gca,'yscale','log')
    set(gca,'xticklabels',[])
    axis square
    xlim([0 4])
    hold on 
end

subplot(2,2,3)
for i=1:size(I2,1)
    plot(1:3,D01_dat(i,3:5),'o','linewidth',2,'MarkerSize',6); drawnow
    set(gca,'yscale','log')
    set(gca,'xticklabels',[])
    axis square
    xlim([0 4])
    hold on 
end

subplot(2,2,4)
for i=1:size(I3,1)
    plot(1:3,D1_dat(i,3:5),'o','linewidth',2,'MarkerSize',6); drawnow
    set(gca,'yscale','log')
    set(gca,'xticklabels',[])
    axis square
    xlim([0 4])
    hold on 
end

% 
% %%%%%% D=1. param_small run only for Dvec=1;
% 
% for i=1:10
%     plot(1:3,param_small(i,3:5),'o'); drawnow
%     set(gca,'yscale','log')
%     set(gca,'xticklabels',[])
%     axis square
%     xlim([0 4])
%     hold on 
% end


























