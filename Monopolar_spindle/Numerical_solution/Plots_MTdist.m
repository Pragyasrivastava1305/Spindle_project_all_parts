set(0,'DefaultAxesFontSize',14,'DefaultTextFontSize',14);
    
           
        
            %%%%%%%%% CONCENTRATION,NUCLEUS, SPINDLE AND MT
            subplot(2,2,1)
            hold off
  
            plot(x,cL_new,'-','linewidth',1)
            hold on
            dsize=60;
            s1=scatter(xN_per(i),0); hold on; s1.Marker='o';     %%%% nucleus position
            set(s1,'SizeData',dsize,'MarkerFaceColor',[0 0.75 0],'MarkerEdgeColor',[0 0.75 0])

            s2=scatter(xC_per(i),0); hold on; s2.Marker='o';     %%%% centrosome position
            set(s2,'SizeData',dsize,'MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0.8 0 0])
  
  
%             %%%%%% MT extent
%             s3=scatter(x(Ipos),zeros(1,size(Ipos,2))); s3.Marker='.'; %% positive arm
%             set(s3,'SizeData',dsize,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
%     
%             s4=scatter(x(Ineg),zeros(1,size(Ineg,2))); s4.Marker='.'; %% negative arm
%             set(s4,'SizeData',dsize,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
%     
    
            xlabel(['Distance along periphery']); ylabel(['c_L(x)'])
            xlim([0 L]); 
            title([' time=',num2str(time(i))]); 

%             axis square
            drawnow  
  
    
            %%%%%% On and off rates 
            subplot(2,2,2)
            hold off
            s0=scatter(xper,Pdfper); hold on;  s0.Marker ='o';
            set(s0,'SizeData',0.05*dsize,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
  
            s1=scatter(xC_per(i),0) ;   s1.Marker = 'o'; 
            set(s1,'SizeData',dsize,'MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0 0 0])
  
            s2=scatter(xN_per(i),0) ;   s2.Marker = 'o'; 
            set(s2,'SizeData',dsize,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0 0])
  
            box on
  
            xlim([0 x(end)]);    ylim([0 max(max(Pdfper))+0.1]); 
            
 
            xlabel('Distance along cell periphery (\mum)')
            ylabel('Pdf of MT ends')
            title([' time=',num2str(time(i))])
            drawnow
  
            %%%% mean velocity and pulling forces   
            subplot(2,2,3)
            hold off
            plot(time(1:end-1),fact/2/taua,'bo','linewidth',1.5,'MarkerSize',2); hold on
            plot(time(1:end),vmean,'k','linewidth',1.5,'MarkerSize',15); 
            xlim([0 tmax]);  
%             axis square; 
            xlabel(['t(min)']);  
            title([' time=',num2str(time(i))])
            ylabel(['Velocity of mean position'])
            drawnow  
    
    
    
            %%%% mean and relative velocities   
            subplot(2,2,4)
            hold off
            plot(time(1:end-1),(fact/taua)-2*(fspr/tauc),'ro','linewidth',1.5,'MarkerSize',2); hold on
            plot(time(1:end),vrel,'k','linewidth',2,'MarkerSize',15); 
            xlim([0 tmax]);  
%             axis square; 
            ylim([-1,1])
            xlabel(['t(min)']); 
            ylabel(['Relative velocity'])
             title([' time=',num2str(time(i))])
            drawnow  