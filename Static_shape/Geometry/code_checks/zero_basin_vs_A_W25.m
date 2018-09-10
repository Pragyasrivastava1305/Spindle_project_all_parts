close all;
rn=8; rs = 1;
plot_mt = 0; plot_ran= 0; plot_opt=0; movie_opt=0;
plt_init_fin = 0; plt_traj = 1; plt_scat =0; 
plt_hist =0; 
nbin = 18; 
phi_initial = round(linspace(0,pi/2,45),rn); 
Asp_vec = 1.2:0.2:5;
b_vec = 12.5; 

tmax=100; dtreg = 0.5; 
t_reg=0:dtreg:tmax; 
Lm = 15; len_sp = 7 ; mu0 = 0.00063;  
ArmR = 9; RanR = 12; koffout = 0.02; koffin =0.1; 

zero_basin = zeros(size(Asp_vec,2),1);
zero_basin(1) = 1;
ang_init = rad2deg(phi_initial);
Npoint =45; iter =1;

dir = 'W_25_sim_data'; 
hist_map = zeros(size(Asp_vec,2),nbin); 
b = 12.5

for iAsp=1:size(Asp_vec,2)
    Asp = Asp_vec(iAsp); 
    a = b*Asp_vec(iAsp);                %%%% semi-major and minor axes  
    shape = @(th) a./sqrt(a^2/b^2 - (a^2/b^2-1)*cos(th).^2);
    th = round(linspace(0,2*pi,360*25),rn);
    dth = round(th(2)-th(1),3);

    rth = shape(th);
    drdth = -a*((a/b)^2 -1)*cos(th).*sin(th)./(((a/b)^2 - ((a/b)^2-1)*cos(th).^2).^1.5);
    fun = sqrt(rth.^2 + drdth.^2); 
    s(1) = 0;

    for ith =2:size(th,2)
         s(ith) = round(trapz(th(1:ith),fun(1:ith)),1); 
    end

    xunit = shape(th).* cos(th) ;
    yunit = shape(th).* sin(th) ;

    %%%%%% system and coords %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% rounding place of L and dx should be same. 
    dx=0.1; 
    L = round(s(end),rs);  
    x = 0:dx:L; 
    N = size(x,2);
    c_unif = 30*ones(1,size(x,2)); 
    phi_sp = load(fullfile(dir,['Npoint=',num2str(Npoint),'Asp=',num2str(Asp_vec(iAsp))...
                                    ,'_Lm=',num2str(Lm),'Phi_dyn_a',num2str(a),'.dat'])); 
    phi_nuc = mod(phi_sp+pi/2,2*pi); 
    figure(iAsp)
    for it =1:Npoint
        omega_init(it,iAsp) =  (phi_sp(it,2)-phi_sp(it,1))/dtreg;
        if plt_traj ==1 
        plot(t_reg,rad2deg(phi_sp(it,:)),'linewidth',1.5); hold on; drawnow
        end
    end
    
    grid on; 
    box on 
    axis square
    set(gca, 'Fontsize',14)
    xlabel('Time','Fontsize',16)
    ylabel('Spindle angle', 'Fontsize',16)
    title(['Aspect Ratio=', num2str(Asp),' width =',num2str(2*b)])
    
    zero_basin(iAsp) = 1+size(find(omega_init(:,iAsp) < 0),1); 
    phi_bound(iAsp) = ang_init(zero_basin(iAsp));
 
    [ycnt,xbin] = hist(rad2deg(phi_sp(:,end)),nbin);
    hist_map(iAsp,:) =  ycnt; 
end
 
if plt_scat ==1
for iAsp = 1:size(Asp_vec,2)
    Ilong =  find(ang_init <= phi_bound(iAsp));
    Ishort= find(ang_init > phi_bound(iAsp));

    scatter(Asp_vec(iAsp)*ones(size(Ilong,2),1),ang_init(Ilong),'MarkerFaceColor','r','MarkerEdgeColor','r',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4); 
    scatter(Asp_vec(iAsp)*ones(size(Ishort,2),1),ang_init(Ishort),'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4); 

hold on 
drawnow
end
end  
 
 
hist_percent = hist_map*100/size(phi_initial,2); 
 if plt_hist ==1
     imagesc(Asp_vec,xbin,hist_percent')
 axis square; axis xy
 set(gca,'Fontsize',15) 
 xlabel('Aspect ratio','Fontsize',18)
 ylabel('Final spindle angle','Fontsize',18)
 
 end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 






