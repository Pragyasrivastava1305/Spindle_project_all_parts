%%%%% Code to generate the phase diagram in the space of geometrical
%%%%% parameters of the system: eccentricity, ecc, half width L, length of
%%%%% the spindle and length of the MTs. lengths are scaled by the half width of the cell.  
tic
%close all
figure;
rn = 8; rs = 1; 
plot_mt = 0; plot_ran= 0; plot_opt=0; movie_opt=0;
Npoint = 90;
opt_disp = 0; 
opt_mov = 0; 
ms = 3; lw = 2; 
xl_el = -25; xh_el= 25;
yl_el = -25; yh_el= 25;
ws  = 0.2;

Lm_vec = 0.1:0.1:2;
ls_vec = 0.1:0.1:1;
% Lm = 15/7 ; 
% len_sp = 1;

%%%  vector of initial spindle orientations 
phi_initial = round(linspace(0,pi/2,Npoint),rn); 
%phi_initial = round(linspace(0.5479,0.5479,Npoint),rn); 
%%%% another specification of initial orientation vector that is finer
%%%% around the unstable FP.

% phi_initial = [round(linspace(0,0.68,20),rn) round(linspace(0.68,0.72,Npoint),rn) round(linspace(0.72,pi/2,20),rn)] 
% mkdir W=25_no_ran
% dir = ['W=',num2str(25),'_no_ran']; 

tmax=50; dtreg = 0.5; 
t_reg=0:dtreg:tmax; 
b = 1;%%%% same as the parameter L of the supplementary
%A = Asp_vec(iAsp); 
ecc = 0.87;
a = b/sqrt(1-ecc^2); 
shape = @(th) b./sqrt(1-(ecc^2)*cos(th).^2);
th = round(linspace(0,2*pi,360*30),rn);
dth = 2*pi/size(th,2);

rth = shape(th);
drdth = -a*((a/b)^2 -1)*cos(th).*sin(th)./(((a/b)^2 - ((a/b)^2-1)*cos(th).^2).^1.5);
fun = sqrt(rth.^2 + drdth.^2); 
s(1) = 0;

for ith =2:size(th,2)
    s(ith) = round(trapz(th(1:ith),fun(1:ith)),rs); 
end

xunit = shape(th).* cos(th) ;
yunit = shape(th).* sin(th) ;

%%%%%% system and coords %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% rounding place of L and dx should be same. 
dx=0.1; 
L = round(s(end),1);  
x = 0:dx:L; 
N = size(x,2);

% iLm =10
% ils =3
for iLm =1: size(Lm_vec,2)    %%%% loop over MT length
for ils = 1: size(ls_vec,2)   %%%%  loop over spindle length 
 [iLm ils]

Lm = Lm_vec(iLm); 
len_sp = ls_vec(ils); 

mu0 = 0.063;  
ArmR = 0.9; RanR = 1.2; koffout = 0.02; koffin =0.1; 

c_unif = ones(1,size(x,2)); 
phi_sp = zeros(size(phi_initial,2),size(t_reg,2)); 
phi_nuc = zeros(size(phi_initial,2),size(t_reg,2));
torq_dens = zeros(size(s,2),size(t_reg,2)-1);
iphi = 5; 

for iphi = 1:size(phi_initial,2)     %%%%  Iteration on all initial angles. 
  %  iphi
phi_sp(iphi,1) = phi_initial(iphi);
phi_nuc(iphi,1) = phi_sp(1) + pi/2;
tic 

%%%%%  initialise movie if opt_mov ==1
if opt_mov ==1 
    mov = VideoWriter(['Movie_phi_in=',num2str(phi_sp(1)),'lc=',num2str(len_sp),...
                        'Lm=',num2str(Lm),'b=',num2str(b),'e=',num2str(ecc),'.avi']); 
end

% phi_sp(1)
for iter = 1:size(t_reg,2)-1       %%%% Iteration over time
%     if iter ==2
%         break
%     end 
%     iter
    %[iphi, iter]
    [torque(iphi,iter),tdens1,tdens2,MT_ang1,MT_ang2] = Torque_static_pinned_Aug18(a,b,th,s,dth,dx,rn,rs...
                                                            ,phi_sp(iphi,iter),len_sp,Lm,x,L,c_unif);
%     ang_chk1(iphi) = MT_ang1(1);
%     ang_chk2(iphi) = MT_ang1(2);
%     size_reg(iphi) = size(MT_ang1,2);
    %%%% new evol
    phi_sp(iphi,iter+1) = round(phi_sp(iphi,iter) + mu0*torque(iphi,iter)*dtreg,rn);
    phi_nuc(iphi,iter+1) = phi_sp(iphi,iter+1) + pi/2;
         
    if opt_disp ==1 
        hold off
        plot(xunit,yunit,'k','linewidth',1.5)
        axis equal
        hold on;
        %%% position and orientation of nucleus
        if (cos(phi_nuc(iphi,iter)) > 0 || phi_nuc(iphi,iter) ==pi/2)
            ang1 = mod(phi_nuc(iphi,iter),2*pi); 
        elseif (cos(phi_nuc(iphi,iter)) < 0 || phi_nuc(iphi,iter) ==3*pi/2)
            ang1 = mod(phi_nuc(iphi,iter)+pi,2*pi); 
        end

        N1x = ArmR*cos(ang1); N1y = ArmR*sin(ang1);
        N2x = ArmR*cos(ang1+pi); N2y = ArmR*sin(ang1+pi);

        %%%% position of centrosomes and orientation of spindle. 
        if (cos(phi_sp(iphi,iter))>0 || phi_sp(iphi,iter)== pi/2)
            phi1 = phi_sp(iphi,iter); phi2 = mod(phi_sp(iphi,iter) + pi,2*pi); 
        elseif (cos(phi_sp(iphi,iter))<0 || phi_sp(iphi,iter) == 3*pi/2)
            phi1 = mod(phi_sp(iphi,iter) + pi,2*pi); phi2 = mod(phi_sp(iphi,iter),2*pi);
        end

        C1x = len_sp*cos(phi1); C1y = len_sp*sin(phi1);
        C2x = len_sp*cos(phi2); C2y = len_sp*sin(phi2);

        ang_s = phi1+pi/2; 
        xTriang = [C1x (ArmR-0.2)*cos(ang_s) C2x (ArmR-0.2)*cos(ang_s+pi)];
        yTriang = [C1y (ArmR-0.2)*sin(ang_s) C2y (ArmR-0.2)*sin(ang_s+pi)]; 
        fill(xTriang,yTriang,'b'); alpha(0.1)

        hold on

        %%%% plot nucleus and its orientation as an ellipse; 
        Da= ArmR; Db = ws; 
        shape = @(th) Da./sqrt(Da^2/Db/Db - (Da^2/Db/Db-1)*cos(th-ang1).^2);
        xdna = shape(th).*cos(th); 
        ydna = shape(th).*sin(th);
        fill(xdna,ydna,'g'); %alpha(0.25)
        circle(C1x,C1y,Lm,lw,'b');   %%% ran discs 
        circle(C2x,C2y,Lm,lw,'b');  
        
        % %%%%% plot astral MT regions on periphery
        if (size(MT_ang1,2)<1)
            %'no torque'
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
            %hold on; plot(xunit(mt_reg2),yunit(mt_reg2),'bo','Markersize',3);
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
           % hold on; plot(xunit(mt_reg3),yunit(mt_reg3),'bo',...
            %                xunit(mt_reg4),yunit(mt_reg4),'bo','Markersize',3)
        end        
        
        title(['frame=',num2str(iter)])
        drawnow
        %%%%  add frame to move if opt_mov ==1  
        if opt_mov ==1
            %%%% add frame to movie
            F = getframe(gcf);
            open(mov);
            writeVideo(mov,F);
        end
    end
  
end
if opt_mov ==1
close(mov)
end
run_one = toc;
% figure(1)
% plot(phi_sp(iphi,:),'.','linewidth',1.5); 
% hold on
% drawnow;

end
close(gcf)
%%%%  save the workspace :  these are not heavy simulations so it is fine
%%%%  to save whole workspace. 
save(['b=',num2str(b),'ecc=',num2str(ecc),'lc=',num2str(len_sp),'Lm=',num2str(Lm),'.mat']);

end     %%%% end of loop over  lc
end     %%%% end of loop over  lm 




