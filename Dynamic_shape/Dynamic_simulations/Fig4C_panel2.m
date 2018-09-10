%%% This code produces Fig. 4D. Some manual edits are needed on the inset
%%% figure. 
%%%% load workspace
load Dyn_sim_Phi_in=70_mu=0.02_dx1=0.01.mat	
%Dyn_sim_Phi_in=70_mu=0.007_dx1=0.002.mat

tpnts = [1 13 19 61 201]; 

for ifr =1:size(tpnts,2)
    plot(xgrid(tpnts(ifr),:),cL_arr(:,tpnts(ifr)),'linewidth',2); hold on
    set(gca,'fontsize',16)
    xlim([0 L(1)]);  ylim([0 6]);
    xlabel('s (\mum)','fontsize',18)
    ylabel('\bar{c}(s)','fontsize',18)
    legend('t=0(NEP)','t=6 min','t=9 min','t=30 min','t=100 min','Location','northwest')
    %axis square
end
%%%% plot spindle angle in time and resposition it manually
axes('Position',[.7 .7 .2 .2])
plot(tarray,rad2deg(phi_sp_reg),'k','linewidth',1.5)
set(gca,'fontsize',12)
xlabel('Time (min)','fontsize',14)
ylabel('Spindle angle (degrees)','fontsize',14)
box on



