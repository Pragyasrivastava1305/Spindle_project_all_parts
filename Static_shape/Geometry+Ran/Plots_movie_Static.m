% close all; clear all
% rn=8; rs = 1; Taud =10;
% plot_mt = 0; plot_ran= 0; plot_opt=0; movie_opt=0;
% plot_phi_dyn =1 ;
% 
ms = 3; lw = 2; 
 xl_el = -25; xh_el= 25;
 yl_el = -25; yh_el= 25;
 ws  = 2;
% 
% 
% phi_initial = round(linspace(0,pi/2,45),rn); 
% Asp_vec = 1.1:0.1:2;
% b_vec = [10 11.25 12.5]; 
% 
% tmax=60; dtreg = 1; 
% t_reg=0:dtreg:tmax; 
% Lm = 15; len_sp = 7 ; mu0 = 0.00063;  
% ArmR = 9; RanR = 12; koffout = 0.02; koffin =0.1; 
% 
% zero_basin = zeros(size(Asp_vec,2),1);
% zero_basin(1) = 1;
% ang_init = rad2deg(phi_initial);
% Npoint =45; iter =1;
% phi_sp = zeros(size(phi_initial,2),size(t_reg,2)); 
% 
% basin_mat = zeros(size(Asp_vec,2),size(phi_initial,2));
% 
% iw = 1; iAsp =9; iphi = 33; 
% b = b_vec(iw); 
% Asp = Asp_vec(iAsp); 
% a = b*Asp_vec(iAsp);   
% 
% %%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% full_data=load(['iPhi=',num2str(iphi),'RanR=',num2str(RanR),'_Taud=',num2str(Taud)...
%      ,'_A=',num2str(Asp),'_b=',num2str(b),'_Lm=',num2str(Lm),'.csv']);
% 
% cL_arr = full_data(1:end-1,:);
% phi_sp_reg =  full_data(end,:);
% mov = VideoWriter(['With_Ran_b=',num2str(b),'_Asp=',num2str(Asp_vec(iAsp)),'_iphi=',num2str(iphi),'.avi']);
% 
% 
% 
% %%%%%%%%%%%%%%%% specify shape, length, space grids. 
%     shape = @(th) a./sqrt(a^2/b^2 - (a^2/b^2-1)*cos(th).^2);
%     th = round(linspace(0,2*pi,360*20),rn);
%     dth = round(th(2)-th(1),3);
% 
%     rth = shape(th);
%     drdth = -a*((a/b)^2 -1)*cos(th).*sin(th)./(((a/b)^2 - ((a/b)^2-1)*cos(th).^2).^1.5);
%     fun = sqrt(rth.^2 + drdth.^2); 
%     s(1) = 0;
% 
%     for ith =2:size(th,2)
%          s(ith) = round(trapz(th(1:ith),fun(1:ith)),1); 
%     end
% 
%     xunit = shape(th).* cos(th) ;
%     yunit = shape(th).* sin(th) ;
% 
%     %%%%%% system and coords %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%% rounding place of L and dx should be same. 
%     dx=0.5; 
%     L = round(s(end),rs);  
%     N = size(full_data,1)-1;
%     x = linspace(0,L,N); 

xinterp = linspace(0,L,size(xunit,2)); 
%phi_nuc_reg= mod(phi_sp_reg+pi/2,2*pi);

%for iframe = 1:size(t_reg,2)
% if iframe >1
%     break 
% end 
% cL_new = cL_arr(:,iframe); 
    
%%% position and orientation of nucleus
if (cos(phi_nuc_reg(iframe)) > 0 | phi_nuc_reg(iframe) ==pi/2)
    ang1 = mod(phi_nuc_reg(iframe),2*pi); 
elseif (cos(phi_nuc_reg(iframe)) < 0 | phi_nuc_reg(iframe) ==3*pi/2)
    ang1 = mod(phi_nuc_reg(iframe)+pi,2*pi); 
end

N1x = ArmR*cos(ang1); N1y = ArmR*sin(ang1);
N2x = ArmR*cos(ang1+pi); N2y = ArmR*sin(ang1+pi);

%%%% position of centrosomes and orientation of spindle. 
if (cos(phi_sp_reg(iframe))>0 | phi_sp_reg(iframe)== pi/2)
    phi1 = phi_sp_reg(iframe); phi2 = mod(phi_sp_reg(iframe) + pi,2*pi); 
elseif (cos(phi_sp_reg(iframe))<0 | phi_sp_save(iframe) == 3*pi/2)
    phi1 = mod(phi_sp_reg(iframe) + pi,2*pi); phi2 = mod(phi_sp_reg(iframe),2*pi);
end

C1x = len_sp*cos(phi1); C1y = len_sp*sin(phi1);
C2x = len_sp*cos(phi2); C2y = len_sp*sin(phi2);


% figure(1)
% subplot(1,2,1)

%%%%% plot ellipse 
% plot(xunit,yunit,'ro','Markersize',ms);

ang_s = phi1+pi/2; 
xTriang = [C1x (ArmR-2)*cos(ang_s) C2x (ArmR-2)*cos(ang_s+pi)];
yTriang = [C1y (ArmR-2)*sin(ang_s) C2y (ArmR-2)*sin(ang_s+pi)]; 
fill(xTriang,yTriang,'b'); alpha(0.1)

hold on

%%%% plot nucleus and its orientation as an ellipse; 
Da= ArmR; Db = ws; 
shape = @(th) Da./sqrt(Da^2/Db/Db - (Da^2/Db/Db-1)*cos(th-ang1).^2);
xdna = shape(th).*cos(th); 
ydna = shape(th).*sin(th);
fill(xdna,ydna,'g'); %alpha(0.25)

hold on;
 
%%%% subplot 2 :  concentration of LGN on contour 
% subplot(1,2,2)
cinterp = interp1(x,cL_arr(:,iframe),xinterp); 
%%% define scaling function to assign colors 
if max(cinterp) == min(cinterp)
    y = cinterp/max(cinterp) ;
else
    y = (cinterp-min(cinterp))/(max(cinterp)-min(cinterp)) ;
end

for ic = 1:size(xunit,2)
    plot(xunit(ic),yunit(ic),'o','MarkerFaceColor',[y(ic) 0 0 ],'MarkerEdgeColor',[y(ic) 0 0],'MarkerSize',ms); 
%     hold on
    
end

axis equal; 
xlim([xl_el xh_el]);  ylim([yl_el yh_el])
box on
set(gca,'Fontsize',13)
title(['time=',num2str(t_reg(iframe))])
xlabel('x','Fontsize',16); 
ylabel('y','Fontsize',16)

hold off

% %%%% add frame to movie
% F = getframe(gcf);
% open(mov);
% writeVideo(mov,F);
% 
% 
% 
% %end
% 
% close(mov)




