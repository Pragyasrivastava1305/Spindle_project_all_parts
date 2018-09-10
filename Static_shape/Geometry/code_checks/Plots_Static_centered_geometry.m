%%%% makes plots. Not defined as a function because there can be arbitrary
%%%% number of changes to be made depending on the final settings required.

ms = 3; lw = 2; dsize = 60;
xl_el = -30; xh_el= 30;
yl_el = -20; yh_el= 20;
col = [0.8 0.8 0.8];



% for i = 1:length(x)
%     hold on
%     plot(xunit,yunit,'','Markersize','.','Color',[1  abs(y(i)) 0])
% end



%%%% position and orientation of nucleus 
if (cos(phi_nuc(iphi,iter)) > 0 )
    ang1 = mod(phi_nuc(iphi,iter),2*pi); 
elseif (cos(phi_nuc(iphi,iter)) < 0 )
    ang1 = mod(phi_nuc(iphi,iter)+pi,2*pi); 
end
% 
% if ( mod(phi_nuc(iphi,iter),pi/2) <= dth )
%     ang1 = pi/2;
% end

N1x = ArmR*cos(ang1); N1y = ArmR*sin(ang1);
N2x = ArmR*cos(ang1+pi); N2y = ArmR*sin(ang1+pi);

%%%% position of centrosomes and orientation of spindle. 
if (cos(phi_sp(iphi,iter))>0)
    phi1 = phi_sp(iphi,iter); phi2 = mod(phi_sp(iphi,iter) + pi,2*pi); 
elseif (cos(phi_sp(iphi,iter))<0)
    phi1 = mod(phi_sp(iphi,iter) + pi,2*pi); phi2 = mod(phi_sp(iphi,iter),2*pi);
end

if ( mod(phi_sp(iphi,iter),pi/2) <=  dth && phi_sp(iphi,iter)~=0)
    phi1 = pi/2; phi2 = 3*pi/2; 
end


C1x = len_sp*cos(phi1); C1y = len_sp*sin(phi1);
C2x = len_sp*cos(phi2); C2y = len_sp*sin(phi2);


% subplot(1,2,1)

%%%%% plot ellipse 
plot(xunit,yunit,'ro','Markersize',ms);

%%%% plot nucleus and its orientation
hold on; plot([N1x,N2x],[N1y,N2y],'og-','linewidth',lw,'Markersize',3*ms...,
    ,'Markeredgecolor','k','Markerfacecolor','g');


if plot_ran ==1
hold on; circle(N1x,N1y,RanR,lw,col); circle(N2x,N2y,RanR,lw,col);

%%%% plot inhibited regions as black 
if size(inhib1,2) <1
    
    
elseif (size(inhib1,2)==2)
    reg1 = (th >= inhib1(1) & th < inhib1(2));  %%% inhib
    reg2 = (th >= inhib2(1) & th < inhib2(2));  %%% inhib
    if (inhib1(1) > inhib1(2))
        reg1 = ( th <= inhib1(2) | th >= inhib1(1));
    end
    hold on; 
    plot(xunit(reg1),yunit(reg1),'ko',xunit(reg2),yunit(reg2),'ko','Markersize',ms); 
    
 elseif (size(inhib1,2)>2)  
    reg1 = (th >= inhib1(1) & th < inhib1(2));  %%% inhib
    reg2 = (th >= inhib1(3) & th < inhib1(4));  %%% inhib
    reg3 = (th >= inhib2(1) & th < inhib2(2));  %%% inhib
    reg4 = (th >= inhib2(3) & th < inhib2(4));  %%% inhib
    
    hold on; 
    plot(xunit(reg1),yunit(reg1),'ko',xunit(reg2),yunit(reg2),'ko',xunit(reg3),yunit(reg3),'ko',...
         xunit(reg4),yunit(reg4),'ko','Markersize',ms);
end
end


%%%% plot centrosomes and spindle
hold on; plot([C1x,C2x],[C1y,C2y],'sq-','Color','[1 0.6 0]','linewidth',lw,'Markersize',4*ms...
         ,'Markeredgecolor','k','Markerfacecolor','[1 0.6 0]'); 
     
if plot_mt ==1

[torque(iter),tdens1,tdens2,MT_ang1,MT_ang2] = Torque_static_pinned_Final_Apr18(a,b,th...
                ,s,dth,rn,rs,phi_sp(iphi,iter),len_sp,Lm,x,L,c_unif);
    
hold on; circle(C1x,C1y,Lm,lw,col); 
circle(C2x,C2y,Lm,lw,col);  

% %%%%% plot astral MT regions 
if (size(MT_ang1,2)<1)
    'no torque'
    
elseif (size(MT_ang1,2)==2)
    
    mt_reg1 = (th >= MT_ang1(1) & th < MT_ang1(2)); 
    mt_reg2 = (th >= MT_ang2(1) & th < MT_ang2(2));  
   
    if (MT_ang1(1) > MT_ang1(2))
    mt_reg1 = ( th <= MT_ang1(2) | th >= MT_ang1(1));
    end   
    
    if (MT_ang2(1) > MT_ang2(2))
    mt_reg2 = ( th <= MT_ang2(2) | th >= MT_ang2(1));   
    end
    %%% overlap region
    hold on; plot(xunit(mt_reg1),yunit(mt_reg1),'bo','Markersize',3);
    hold on; plot(xunit(mt_reg2),yunit(mt_reg2),'bo','Markersize',3);
 elseif (size(MT_ang1,2)>2)
    
    mt_reg1 = (th >= MT_ang1(1) & th < MT_ang1(2));  
    mt_reg2 = (th >= MT_ang1(3) & th < MT_ang1(4));  
    

    if (sin(MT_ang1(1))*sin(MT_ang1(2))) < 0
        mt_reg1 = (th >= MT_ang1(1) | th < MT_ang1(2));
    end
    
    
    if (sin(MT_ang1(3))*sin(MT_ang1(4))) < 0
        mt_reg2 = (th >= MT_ang1(3) | th < MT_ang1(4));
    end    
    
    mt_reg3 = (th >= MT_ang2(1) & th < MT_ang2(2));  
    mt_reg4 = (th >= MT_ang2(3) & th < MT_ang2(4));  
    
    hold on; plot(xunit(mt_reg1),yunit(mt_reg1),'bo',xunit(mt_reg2),yunit(mt_reg2)...
                    ,'bo','Markersize',3)
    hold on; plot(xunit(mt_reg3),yunit(mt_reg3),'bo',...
                    xunit(mt_reg4),yunit(mt_reg4),'bo','Markersize',3)
                
end

end
axis equal; 
xlim([xl_el xh_el]);  ylim([yl_el yh_el])
xlabel('x','Fontsize',14); 
ylabel('y','Fontsize',14); 
% title(['time=',num2str(t_reg(iter))])
hold off
%%%% plot concentration



