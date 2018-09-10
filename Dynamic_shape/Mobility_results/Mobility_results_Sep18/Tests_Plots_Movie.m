opt_movie =1;   iphi =9;
opt_plt =0; Tround =10; dx =0.02;
T_total=450;  cl_lo = 5;  cl_hi = 40; 

dir = 'Ttotal=450';
phi_vec  = deg2rad(linspace(0,90,10)); 
if opt_plt ==1
for iphi =1:size(phi_vec,2)
   load(fullfile(dir,['Phi=',num2str(phi_vec(iphi)),'_Tround=',num2str(Tround),'dx=',num2str(dx(1)),'.mat'])); 
   %load(); 
   figure(1)
   plot(tarray,cell_len,'linewidth',1.5); hold on; drawnow
   plot(tarray,cell_wid,'linewidth',1.5)
   figure(2)
   plot(tarray,rad2deg(phi_sp_reg),'linewidth',1.5); hold on; drawnow
   phi_fin(iphi) = rad2deg(phi_sp_reg(end)); 
   phi_20(iphi) = rad2deg(phi_sp_reg(41)); 
end

end


if opt_movie ==1 
load(fullfile(dir,['Phi=',num2str(phi_vec(iphi)),'_Tround=',num2str(Tround),'dx=',num2str(dx(1)),'.mat'])); 
         
mov = VideoWriter(['Dynamic_simulation_phi=80_dx=02_c0=30_circle.avi']);
% mov = VideoWriter(['Conc_diff_phi=80_dx=02_c0=30.avi']);

for iframe = 1:size(tarray,2)-1

    
shape = @(th) cell_a(iter)./sqrt((cell_len(iframe)/cell_wid(iframe))^2 ...
                - ((cell_len(iframe)/cell_wid(iframe))^2-1)*cos(th).^2);
th = round(linspace(0,2*pi,360*10),rn);
dth = round(th(2)-th(1),3);
rth = shape(th);
den = ((cell_len(iframe)/cell_wid(iframe))^2- ((cell_len(iframe)/cell_wid(iframe))^2-1)*cos(th).^2).^1.5;
drdth = -cell_len(iframe)*((cell_len(iframe)/cell_wid(iframe))^2 -1)*cos(th).*sin(th)./den;
fun = sqrt(rth.^2 + drdth.^2); 
s_in(1) = 0;

for ith =2:size(th,2)
      s_in(ith) = round(trapz(th(1:ith),fun(1:ith)),rs); 
end

tpnt = iframe; 
conc = cL_arr(:,iframe); 
tim = tarray(iframe); 

if tarray(iframe) < 10
    Dynamic_plots_running_sim;
else
    Dynamic_plots_running_sim_with_angle;
end    
    
% if tarray(iframe) > 10
%     plot(cL_arr(:,iframe+1)-cL_arr(:,iframe),'linewidth',1.5)
%     title(['time=',num2str(tarray(iframe))])
%     ylim([-1 1])
% end

 F = getframe(gcf);
 open(mov); 
 writeVideo(mov,F);    
 end


close(mov)

end














