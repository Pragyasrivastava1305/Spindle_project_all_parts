set(gcf,'renderer','painters');
% clear all
% close all; 
ms = 4; lw = 2; 
ws  = 0.5; 
plt_movie =0; 
load Dyn_sim_Phi_in=70_mu=0.02_dx1=0.01.mat	

cl_lo = 0; cl_hi = 6;

%%%% specify time point for the tpnt. 
tpnt =1;
xl_el = -22; xh_el= 22;
yl_el = -15; yh_el= 15;
    
conc = cL_arr(:,tpnt); 
tim =  tarray(tpnt); 

%%% define ellipse shape 
ecc = sqrt(1-(cell_wid(tpnt)^2)/cell_len(tpnt)/cell_len(tpnt));
rplot = @(th) cell_wid(tpnt)./sqrt(1-(ecc^2)*cos(th).^2);
xunit= rplot(th).* cos(th) ;
yunit = rplot(th).* sin(th) ;
rth = rplot(th);
drdth = -cell_len(tpnt)*(ecc^2)*cos(th).*sin(th)./((1-(ecc^2)*(cos(th)).^2).^1.5);
fun = sqrt(rth.^2 + drdth.^2); 
s_plot(1) = 0;

for ith =2:size(th,2)
    s_plot(ith) = round(trapz(th(1:ith),fun(1:ith)),rs); 
end

%%%%% get concentration on s. 
xinterp = linspace(0,xgrid(tpnt,end),size(xunit,2)); 
cinterp = interp1(xgrid(tpnt,:),conc,xinterp); 
%%% define scaling function to assign colors 
if max(cinterp) == min(cinterp)
    y = cinterp/max(cinterp) ;
else
    y = (cinterp-min(cinterp))/(max(cinterp)-min(cinterp)) ;
end

    
for ic = 1:size(xunit,2)
    if tpnt ==1
        plot(xunit(ic),yunit(ic),'.','MarkerFaceColor',[0 0 0 ],'MarkerEdgeColor',[0 0 0],'MarkerSize',ms); 
        hold on
    else
        plot(xunit(ic),yunit(ic),'.','MarkerFaceColor',[y(ic) 0 0 ],'MarkerEdgeColor',[y(ic) 0 0],'MarkerSize',ms); 
        hold on
    end
end
   
hold on 

if tim < Tround
        
%%%%%%%%%%%%%% Plot inhibited regions and the changing shape
%         if dna_R < cell_len(tpnt)    
%         plot(xunit,yunit,'r.','Markersize',ms);
%         elseif dna_R >= cell_len(tpnt)
%         plot(xunit,yunit,'g.','Markersize',ms);
%         end
xcircle = (dna_R-4)*cos(th); 
ycircle = (dna_R-4)*sin(th);
hold on; 
%circle(0,0,dna_R-4,lw,'g');
fill(xcircle,ycircle,'b'); alpha(0.65)
hold on;

distance= sqrt((rplot(th).*cos(th)).^2+(rplot(th).*sin(th)).^2);
fun_dist = dna_R - distance; 
%%%%% plot regions for which dist_fun >0; 
% plot(xunit(fun_dist>=0),yunit(fun_dist>=0),'g.','markersize',ms)
%%%% get off and on rates
[suniquep,iuniquep,ichk] = unique(s_plot);
offrate_on_s_p = koffout*ones(size(s_plot,2),1); 
offrate_on_s_p(dist_fun >=0) = koffin; 
koff_plot = interp1(suniquep,offrate_on_s_p(iuniquep),x); 

 axis equal
 xlim([xl_el xh_el]); ylim([yl_el yh_el])
 title(['time=',num2str(tim)])
 xlabel('x','Fontsize',16); 
 ylabel('y','Fontsize',16)
 hold off; axis equal
 drawnow
    

else
        
%         xunit= rplot(th).* cos(th) ;
%         yunit = rplot(th).* sin(th) ;
%         plot(xunit,yunit,'r.','Markersize',ms);

        hold all;

        %%%% position of centrosomes and orientation of spindle. 
        if (cos(phi_sp_reg(tpnt))>0 || phi_sp_reg(tpnt)== pi/2)
            phi1 = phi_sp_reg(tpnt); phi2 = mod(phi_sp_reg(tpnt) + pi,2*pi); 
        elseif (cos(phi_sp_reg(tpnt))<0 || phi_sp_save(tpnt) == 3*pi/2)
            phi1 = mod(phi_sp_reg(tpnt) + pi,2*pi); phi2 = mod(phi_sp_reg(tpnt),2*pi);
        end

        C1x = len_sp*cos(phi1); C1y = len_sp*sin(phi1);
        C2x = len_sp*cos(phi2); C2y = len_sp*sin(phi2);

        %%%%%% plot spindle as two triangles 
        ang_s = phi1+pi/2; 
        xTriang = [C1x (ArmR-2)*cos(ang_s) C2x (ArmR-2)*cos(ang_s+pi)];
        yTriang = [C1y (ArmR-2)*sin(ang_s) C2y (ArmR-2)*sin(ang_s+pi)]; 
        fill(xTriang,yTriang,'g'); alpha(0.4)

        %%% position and orientation of nucleus
        if (cos(phi_nuc_reg(tpnt)) > 0 || phi_nuc_reg(tpnt) ==pi/2)
            ang1 = mod(phi_nuc_reg(tnpt),2*pi); 
        elseif (cos(phi_nuc_reg(tpnt)) < 0 || phi_nuc_reg(tpnt) ==3*pi/2)
            ang1 = mod(phi_nuc_reg(tpnt)+pi,2*pi); 
        end

        N1x = ArmR*cos(ang1); N1y = ArmR*sin(ang1);
        N2x = ArmR*cos(ang1+pi); N2y = ArmR*sin(ang1+pi);

        %%%% plot nucleus and its orientation as an ellipse; 
        Da= ArmR; Db = ws; 
        shape_dna = @(th) Da./sqrt(Da^2/Db/Db - (Da^2/Db/Db-1)*cos(th-ang1).^2);
        xdna = shape_dna(th).*cos(th); 
        ydna = shape_dna(th).*sin(th);
        fill(xdna,ydna,'b'); alpha(0.65)


        %%%%%%%% plot inhibited regions and the regions of torque generation
        %%%%%%%% write these as functions later on. 
        %%%%%% plot astral MT regions 
        fun1 = @(th) Lm-sqrt(len_sp^2 - 2*len_sp*rplot(th).*cos(th-phi1)+ rplot(th).^2);
        fun2 = @(th) Lm-sqrt(len_sp^2 - 2*len_sp*rplot(th).*cos(th-phi2)+ rplot(th).^2);
        tvec1 = fun1(th)>=0; tvec2 = fun2(th) >=0; 

%         plot(xunit(tvec1),yunit(tvec1),'b.',...
%              xunit(tvec2),yunit(tvec2),'b.','Markersize',ms); hold on; 

        %%%%%% plot inhibited regions on the periphery
        distance = closest_dist(rplot,th,[N1x N1y],[N2x,N2y]);
        fun_dist = RanR - distance; 
        %plot(xunit(fun_dist>=0),yunit(fun_dist>=0),'g.','Markersize',ms)
        %%%%  calculate koff_plot to plot later.
        dist_PtoN_p = closest_dist(shape,th,[N1x N1y],[N2x,N2y]);
        dist_fun_p = RanR - dist_PtoN_p; 

        %%%% get off and on rates
        [suniquep,iuniquep,ichk] = unique(s_plot);
        offrate_on_s_p = koffout*ones(size(s_plot,2),1); 
        offrate_on_s_p(dist_fun_p >=0) = koffin; 
        koff_plot = interp1(suniquep,offrate_on_s_p(iuniquep),x); 

        axis equal
        xlim([xl_el xh_el]); ylim([yl_el yh_el])
        title(['time=',num2str(tim)])
        set(gca,'fontsize',16)
        xlabel('x','Fontsize',18); 
        ylabel('y','Fontsize',18)
        hold off
        drawnow
        axis equal
        hold off;
      
end
           
saveas(gcf,['T=',num2str(tim),'_snapshot.eps'])
