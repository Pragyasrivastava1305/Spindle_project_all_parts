%%%%% CODE TO GENERATE WIDTH AND AR DEPENDENCE OF THE GEOMETRY DEPENDENT
%%%%% CORRECTION.
%%%%% Simplified code to get only geometry dependent dynamics and steady
%%%%% state of the spindle angle. This is same as the main code but much
%%%%% simplified because reaction diffusion is not implemented.
%%%%% Concentration is specified to be uniform at all times. 
% clear all; 
% close all;  
% dir = 'Only_geometric_static_pinned';
% dir = 'Dependence_on_dphi';
tic
%close all
figure;
rn = 8; rs = 1; 
plot_mt = 0; plot_ran= 0; plot_opt=0; movie_opt=0;
Npoint = 50;
% phi_initial = round(linspace(0,pi/2,Npoint),rn); 
phi_initial = [round(linspace(0,0.68,20),rn) round(linspace(0.68,0.72,Npoint),rn) round(linspace(0.72,pi/2,20),rn)] 
% mkdir W=25_no_ran
% dir = ['W=',num2str(25),'_no_ran']; 

%%%% parameter range for histograms
% wid_vec = 12.5; 
% Asp_vec = [1.2:0.2:5]; 

        
tmax=5; dtreg = 0.05; 
t_reg=0:dtreg:tmax; 
    
a = 20; b =10;
shape = @(th) a./sqrt(a^2/b^2 - (a^2/b^2-1)*cos(th).^2);
th = round(linspace(0,2*pi,360*25),rn);
% dth = round(th(2)-th(1),3);
% dth =  mean(diff(th));
dth = 2*pi/size(th,2);

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
phi_sp = zeros(1,size(t_reg,2)); 
torq_dens = zeros(size(s,2),size(t_reg,2)-1);
torq_densp = zeros(size(s,2),size(t_reg,2)-1);


iphi = 2;
for iphi = 1:size(phi_initial,2)
phi_initial(iphi);
phi_sp(1) = phi_initial(iphi);
phi_spp(1) = phi_initial(iphi);
tic 
% phi_sp(1)
for iter = 1:size(t_reg,2)-1
    %iter 
    %%%% calculate torque
    [torque(iter),tdens1,tdens2,MT_ang1,MT_ang2] = Torque_static_pinned_Final_Apr18(a,b,th,s,dth,rn,rs,phi_sp(iter),len_sp,Lm,x,L,c_unif);
    torq_dens(:,iter) = tdens1+tdens2; 
    
    [torquep(iter),tdensp1,tdensp2,MT_angp1,MT_angp2] = Torque_static_pinned_Dev2(a,b,th,s,dth,rn,rs,phi_spp(iter),len_sp,Lm,x,L,c_unif);
    torq_densp(:,iter) = tdensp1+tdensp2; 
    
%   diff_iter(:,iter) = torq_densp(:,iter)/max(torq_densp(:,iter)) -  torq_dens(:,iter)/max(torq_dens(:,iter)); 
    Tot_old(iter) = trapz(dth*torq_dens(:,iter)); 
    Tot_new(iter) = trapz(dth*torq_densp(:,iter)); 
    
    
%     plot(t)
    %%%% old evol
    phi_sp(iter+1) = round(phi_sp(iter) + mu0*torque(iter)*dtreg,rn);
  
    %if abs(phi_sp(iter+1)- phi_sp(iter)) < dth; phi_sp(iter+1) = phi_sp(iter); end
    
    %%%% new evol
    phi_spp(iter+1) = round(phi_spp(iter) + mu0*torquep(iter)*dtreg,rn);
    %%%% if the change in angle is less than the angular grid spacing, set
    %%%% it to 0. 
    %if abs(phi_spp(iter+1)- phi_spp(iter)) < dth; phi_spp(iter+1) = phi_spp(iter); end
    
end
run_one = toc;
%figure(1)
%plot(phi_sp,'.','linewidth',1.5); 
hold on
plot(t_reg,phi_spp,'.-.','LineWidth',1.5); 
drawnow;

end








