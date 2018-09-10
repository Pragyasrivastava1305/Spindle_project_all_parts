%%%% Code to calculate the noisy dynamics of spindle angle without Ran. LGN
%%%% concentration has been taken to be uniformly high set by the bare kon
%%%% and koff. 
clear all

%%%% Parameters set from monopolar model
rn = 8; rs = 1;
load_sys_params; 

ArmR = params.dna_half_length;
RanR = params.inhibition;
Lm = params.astral_Mts;
len_sp = params.spindle_half_length ; 
D = params.Diffusion; 
mu0 = params.mobility; 
Taud = params.Delay; 
konout = params.kon_far; 
konin = params.kon_near; 
koffout = params.koff_far; 
koffin = params.koff_near; 

%%%%% noise strength and number of iterations
% noi_str_vec = [0.05 0.1 1 5]; 
noi_str_vec = 80; 
noi_iter = 200; 
noi_vec = rand(1,1e8)-0.5;  %%%%% an array containing random numbers with 0 mean and uniform distribution.
phi_initial = 0;

%%%% system specification 
tmax = 500;
a = 12; b = 12; 
shape = @(th) a./sqrt(a^2/b^2 - (a^2/b^2-1)*cos(th).^2);
th = round(linspace(0,2*pi,360*15),rn);
dth = round(th(2)-th(1),3);

rth = shape(th);
drdth = -a*((a/b)^2 -1)*cos(th).*sin(th)./(((a/b)^2 - ((a/b)^2-1)*cos(th).^2).^1.5);
fun = sqrt(rth.^2 + drdth.^2); 
s(1) = 0;

for ith =2:size(th,2)
      s(ith) = round(trapz(th(1:ith),fun(1:ith)),rs); 
end

xunit = shape(th).* cos(th) ;
yunit = shape(th).* sin(th) ;
%%%% plot ellipse as check if needed
% plot(xunit,yunit); axis equal

%%%%% Specify grids
%%%%% Time 
dtc = 0.05; 
dtreg = 1; t_reg = 0:dtreg:tmax;
nsize = size(t_reg,2);
tol_adap = 0.001; 

%%% Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx=0.05;              
L = round(s(end),rs);  
x = 0:dx:L; 
N = size(x,2); 

%%%% define the 0 of torque for given a and b. Its a small discretisation
%%%% error. 
c_unif = 50*ones(1,size(x,2)); 
tol_torque =  abs(Torque_static_pinned_Final_Apr18(a,b,th,s,dth,rn,rs,0,len_sp,Lm,x,L,c_unif));
phi_sp_reg = zeros(noi_iter,nsize); 
phi_sp_reg(:,1) = phi_initial; 

                                
for jnoise = 1:size(noi_str_vec,2)
% jnoise =3;
    noi_str = noi_str_vec(jnoise)
for N_iter = 1: noi_iter
    N_iter 
for T_iter =1:size(t_reg,2)-1
  
[torque(T_iter),tdens1,tdens2,MT_ang1,MT_ang2] = Torque_static_pinned_Final_Apr18(a,b,th,s,dth,...
                                                rn,rs,phi_sp_reg(N_iter,T_iter),len_sp,Lm,x,L,c_unif);
torq_dens(:,T_iter) = tdens1+tdens2; 
    
%%% update angle
eta(T_iter) = noi_str*noi_vec(randperm(length(noi_vec),1));  
phi_sp_reg(N_iter,T_iter+1) = round(phi_sp_reg(N_iter,T_iter) + mu0*torque(T_iter)*dtreg + eta(T_iter), rn) ;
phi_nuc_reg(N_iter,T_iter+1) = round(phi_sp_reg(N_iter,T_iter+1) + pi/2,rn);

end
% plot(t_reg,rad2deg(phi_sp_reg(N_iter,:))) ;  hold on ; drawnow
end  

%%%%% save angle dynamics for all iterations in one file
dlmwrite(['AngleData_No_ran_T=',num2str(tmax),'_Noise=',num2str(noi_str)...
    ,'_a=',num2str(a),'.csv'],phi_sp_reg,'delimiter','\t');    
end


