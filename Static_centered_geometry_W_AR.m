%%%%% CODE TO GENERATE WIDTH AND AR DEPENDENCE OF THE GEOMETRY DEPENDENT
%%%%% CORRECTION.
%%%%% Simplified code to get only geometry dependent dynamics and steady
%%%%% state of the spindle angle. This is same as the main code but much
%%%%% simplified because reaction diffusion is not implemented.
%%%%% Concentration is specified to be uniform at all times. 
clear all; 
close all;  
% dir = 'Only_geometric_static_pinned';
% dir = 'Dependence_on_dphi';
tic
rn = 8; rs = 1; 
plot_mt = 0; plot_ran= 0; plot_opt=0; movie_opt=0;
Npoint = 45;
phi_initial = round(linspace(0,pi/2,Npoint),rn); 
mkdir W=25_no_ran
% dir = ['W=',num2str(25),'_no_ran']; 

%%%% parameter range for histograms
% wid_vec = 12.5; 
% Asp_vec = [1.2:0.2:5]; 

%%%% parameter values to compare with conc and delay checks
Asp_vec = 2 ;
wid_vec = 10; 
for iwid = 1:size(wid_vec,2)
for iAsp = 1:size(Asp_vec,2)
    
    tmax=100; dtreg = 0.5; 
    t_reg=0:dtreg:tmax; 

    %%%%% specify shape of ellipse, calculate perimeter=linear system size. 
    b = wid_vec(iwid);
    a = b*Asp_vec(iAsp);                
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

    Lm = 15; len_sp = 7 ; mu0 = 0.00063;  
    ArmR = 9; RanR = 12; koffout = 0.02; koffin =0.1; 

    c_unif = 50*ones(1,size(x,2)); 
    phi_sp = zeros(size(phi_initial,2),size(t_reg,2)); 
    torq_dens = zeros(size(s,2),size(t_reg,2)-1);

    phi_sp(:,1) = phi_initial; 
    phi_nuc(:,1) = round(phi_sp(:,1) + pi/2,rn);

    for iphi = 1:size(phi_initial,2)
    [iwid iAsp iphi]
    
    if movie_opt ==1
    mov = VideoWriter(['Lm=',num2str(Lm)'...
                    '_a=',num2str(a),'_A=',num2str(A),'_phi=',num2str(round(rad2deg(phi_sp(iphi,1)),1)),'.avi']);
        %  mov=VideoWriter(['Test_phi=',num2str(round(rad2deg(phi_sp(iphi,1)),1)),'.avi']); 
    end

    for iter = 1:size(t_reg,2)-1
    %%%% calculate torque
    [torque(iter),tdens1,tdens2,MT_ang1,MT_ang2] = Torque_static_pinned_Final_Apr18(a,b,th,s,dth,rn,rs,phi_sp(iphi,iter),len_sp,Lm,x,L,c_unif);
    torq_dens(:,iter) = tdens1+tdens2; 
    
    %%% update angle
    phi_sp(iphi,iter+1) = round(phi_sp(iphi,iter) + mu0*torque(iter)*dtreg,rn);
    phi_nuc(iphi,iter+1) = round(phi_sp(iphi,iter+1) + pi/2,rn);
    %%%% Plots and movie
    if plot_opt == 1
        [inhib1,inhib2,koff] = Rates_static_ellipse_Final_Jan18(x,th,s,a,b,phi_nuc(iphi,iter),ArmR,RanR,koffout,koffin); 
        Plots_Static_centered_geometry; drawnow;
        if movie_opt ==1
            F = getframe(gcf);
            open(mov);
            writeVideo(mov,F);
         end
    end
    end
    % % close(mov)
    % dlmwrite(fullfile(dir,['L_m=',num2str(Lm),'Torque_iPhi=',num2str(iphi),'_a=',num2str(a),'b=',num2str(b),'.csv']),torq_dens,'delimiter','\t')   

    end
    
    % plot(t_reg(1:end-1),res_torque(iphi,:),'o-.','linewidth',2,'Markersize',6); 
    % drawnow
    % hold on; 
    dlmwrite(fullfile(dir,['Npoint=',num2str(Npoint),'Asp=',num2str(Asp_vec(iAsp)),'_Lm=',num2str(Lm)...
          ,'Phi_dyn_a',num2str(a),'.mat']),phi_sp,'delimiter','\t')
     
%     dlmwrite(['Npoint=',num2str(Npoint),'Asp=',num2str(Asp_vec(iAsp)),'_Lm=',num2str(Lm)...
%           ,'Phi_dyn_b=',num2str(b),'.dat'],phi_sp,'delimiter','\t') 
%    
    run_shape = toc; 
    

    %%%%%%%  all plots 
    % if plot_opt ==1
    % figure;
    % subplot(1,2,1)
    % for iphi = 1:size(phi_initial,2)
    % plot(t_reg(1:end-1),res_torque(iphi,:),'o-.','linewidth',2,'Markersize',6)
    % hold on; 
    % drawnow
    % xlim([0 t_reg(end)])
    % end
    % 
    % subplot(1,2,2)
    % plot(phi_initial,phi_sp(:,end),'linewidth',2)
    % end 
end 
% clear phi_sp;  clear phi_nuc;

end














