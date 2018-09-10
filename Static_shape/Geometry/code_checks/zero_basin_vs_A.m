close all;
rn=8; rs = 1;
plot_mt = 1; plot_ran= 0; plot_opt=0; movie_opt=0;

phi_initial = round(linspace(0,pi/2,45),rn); 
Asp_vec = 1:0.1:2;
b_vec = [10 11.25 12.5 15]; 

tmax=100; dtreg = 0.1; 
t_reg=0:dtreg:tmax; 
Lm = 15; len_sp = 7 ; mu0 = 0.00063;  
ArmR = 9; RanR = 12; koffout = 0.02; koffin =0.1; 

zero_basin = zeros(size(Asp_vec,2),1);
zero_basin(1) = 1;
ang_init = rad2deg(phi_initial);
Npoint =45; iter =1;

dir = 'W=25_no_ran'; 

 for iw = 1:size(b_vec,2)
% iw = 3; iAsp =1;
for iAsp=1:size(Asp_vec,2)
    
    %%%% specify shape of ellipse, calculate perimeter=linear system size. 
    b = b_vec(iw);  a = b*Asp_vec(iAsp);                %%%% semi-major and minor axes  
    shape = @(th) a./sqrt(a^2/b^2 - (a^2/b^2-1)*cos(th).^2);
    th = round(linspace(0,2*pi,360*20),rn);
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
    
    load(fullfile(dir,['Phi=',num2str(phi_vec(iphi)),'_Tround=',num2str(Tround)...
       '_C=',num2str(c_vec(iclev)),'dx=',num2str(dx1),'.mat'])); 
    phi_sp= load(['Npoint=',num2str(Npoint),'Asp=',num2str(Asp_vec(iAsp)),'_Lm=',num2str(Lm),'Phi_dyn_b=',num2str(b_vec(iw)),'.dat']);
    phi_nuc = mod(phi_sp+pi/2,2*pi); 
%     figure; 
    for it =1:Npoint
        omega_init(it,iAsp) =  (phi_sp(it,2)-phi_sp(it,1))/dtreg;
        plot(t_reg,rad2deg(phi_sp(it,:)),'linewidth',1.5); hold on; drawnow
    end
    
    if iAsp >1
    zero_basin(iAsp) = 1+size(find(omega_init(:,iAsp) < 0),1); 
    figure(iw)
    scatter(Asp_vec(iAsp)*ones(zero_basin(iAsp),1),ang_init(1:zero_basin(iAsp)));
    hold on
    end
    
    phi_bound(iw,iAsp) = ang_init(zero_basin(iAsp));
    
    %%%%  make geometric plots at the first point which converges to pi/2. 
%     figure;
    iphi = zero_basin(iAsp);
%      figure(iw)
%      plot(rad2deg(phi_sp(:,1)),rad2deg(phi_sp(:,end)));  hold on 
     drawnow
    
%      Plots_Static_centered_geometry; axis off;
%      drawnow
%      saveas(gcf,['wid=',num2str(b),'Asp',num2str(iAsp),'.png'])
%     
%    
    
  
% hold on; drawnow   
% end
% figure; 
% plot(Asp_vec,phi_bound(iw,:)); hold on; drawnow; 
% 
end
end


iw =1; 
for iAsp = 2:11;
    Ilong =  find(ang_init <= phi_bound(iw,iAsp));
    Ishort= find(ang_init > phi_bound(iw,iAsp));

    scatter(Asp_vec(iAsp)*ones(size(Ilong,2),1),ang_init(Ilong),'MarkerFaceColor','r','MarkerEdgeColor','r',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4); 
    scatter(Asp_vec(iAsp)*ones(size(Ishort,2),1),ang_init(Ishort),'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.4); 

hold on 
drawnow
end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 






