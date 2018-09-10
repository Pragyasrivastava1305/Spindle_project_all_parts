close all 
clear all

load D=0.01_va=0.9_tauc=0.mat;
%%%%  all the variables are predefined in the mat file, so they do not need
%%%%  to be defined here again. 
mov = VideoWriter(['Mono_Simulation_D=',num2str(0.01),'_dx=',num2str(0.02),...
    '_koffout=',num2str(0.129),'_koffin=',num2str(0.025),'.avi']); 

for ifr = 1: size(t_reg,2)
    hold off
    plot(x,koffin*cL_arr(:,ifr),'linewidth',2,'color',[1 0.5 0.3]); 
    hold on;
    
    dsize=60;
    s1=scatter(xC_per_reg(ifr),0); hold on; s1.Marker='o';  
    set(s1,'SizeData',dsize,'MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0.8 0 0])

    s2=scatter(xN_per_reg(ifr),0); hold on; s2.Marker='o';  
    set(s2,'SizeData',dsize,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0])
    
    
    xlim([0 L]);
    ylim([0 3]); 
    xlabel('x (\mum)','FontSize',16);
    ylabel('LGN profile ','FontSize',16,'Interpreter','tex'); 
    title(['Time=',num2str(t_reg(ifr))],'FontSize',16,'FontWeight','normal');
    legend('c(x,t)/c_{0}^{near} ','Centrosome','DNA')
    set(gca,'FontSize',14)
    
    drawnow
    F = getframe(gcf);
    open(mov);
    writeVideo(mov,F);
end

close(mov)

