%%%% Main testing code for shape and patterning dynamics %%%%%%%%%%%%%%%%%%
%%%% At each iteration before Tround, length of the perimeter changes,
%%%% inhibition range on periphery is calculated with DNA being round and,
%%%% concentration is evolved. At each iteration after Tround, (a) length
%%%% changes, (b) concentration evolves with inhibition range calculated
%%%% with DNA being a line and at a specified angle, (c) torque on the
%%%% spindle is calculated.  In both the cases, at every step,
%%%% concentration needs to be scaled up according to the shrinking factor
%%%% of length.
clear all; 
% mov = VideoWriter(['Dynamic_simulation_test.avi']);
%%%%%% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%%%%%%%%%%%%%%  Fixed PHYSICAL AND NUMERICAL parameters %%%%%%%%%%%%%%%%%%%
rn = 8; rs = 10;

jmu = 6;  %%%% specify jmu and nphi to run the code.
mu0_vec = [1 0.1 0.01 0.02 0.001 1e-5];
mu0 = mu0_vec(jmu);

phi_initial = linspace(0,pi/2,19);
nphi = 2;  

%%%%  to submit on the cluster, 20 points are broken into 4 groups.
%%%%  Depending on nphi, the code selects which group it reads
%%%%  for initial angles. 

%%% Assign physical parameters.
ArmR = 9;    %%% half length of metaphase plate
RanR = 11;   %%% Range of RAN inhibiton
Lm = 15;     %%% Length of MTs
len_sp = 7;  %%% half length of spindle
D = 0.01;    %%% Diffusion constant
%mu0 = 0.02; %%% mobility
Taud = 9;     %%% Delay between NEP and force application on spindle
konout = 1; konin = 1;          %%% On rates outside and inside inhibiton range
koffout = 0.025; koffin = 0.129;  %%% off rates outside and inside. 
tau =1/koffin; 
dna_R = 13;        %%%%%% radius of DNA(9) + inhibition range(4) before t= Tround

T_total = 100;                                    %%%%% Total time
Tround = Taud;                                    %%%%% rounding time
dtreg = 0.5; 
tarray = 0:dtreg: T_total; 
nsize = size(tarray,2); 
dtc = 0.005; tol_adap = 0.001; 

%%%%%% cell_length a*(1-tanh((tnew-C)/tau))+b: obtained from 3 params fit,
%%%%%% fixing b =  22.6250
p1 = 14.27; p2 = 0.1412; p3 = 5.936; p4 = 22.625;
half_cell_len = @(t) 0.5*(p1*(1-tanh((t-p2)/p3))+ p4); 

%%%%% cell aspect ratio
p1 = 0.3965; p2 = 2.706; p3 = 4.866; 
Asp_dyn = @(t) p1*(1-tanh((t-p2)/p3))+1; 

%%%%% cell half length : 3 params fit. b= 22.6250
p1 = 4.326; p2 = -4.9; p3 = 5.46; p4 = 22.6250;
half_cell_wid = @(t) 0.5*(p1*(1-tanh((t-p2)/p3))+ p4); 

%%%%% options to make plots and movies :  
plt_iter = 0;   %%% plot at iterations on adaptive time points
plt_reg =0;     %%% plot at regular time points
plt_movie =0;   %%% movie
ran_inhibition =1;   %%%%  switch Ran inhibition on/off. 
%%%%% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%%%%%%%%%%%%%%%%%%% INITIAL SETUP OF THE SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%
%mov = VideoWriter(['Rounding_with_Ran_movie.avi']);
if nphi ==4; npts = 4; else;  npts =5; end

for iphi =1:npts
[mu0 iphi]
phi_ind = (nphi-1)*5 + iphi;
phi_in = phi_initial(phi_ind);                           %%%%%%%%% ANGLE AT T= TROUND 
iter =1;  jreg =1; dx = 0.02; 
delt(iter) = 0.01;
time(iter) = 0; 

cell_a(iter) = half_cell_len(time(iter));       %%% assign cell half length and width according to experimental fit
cell_b(iter) = half_cell_wid(time(iter));
ecc = sqrt(1-(cell_b(iter)^2)/cell_a(iter)/cell_a(iter));    
 
shape = @(th) cell_b(iter)./sqrt(1-(ecc^2)*cos(th).^2);  %%%%  define ellipse
th = round(linspace(0,2*pi,360*10),rn);
dth = 2*pi/(size(th,2)-1);
rth = shape(th);
drdth = -cell_b(iter)*(ecc^2)*cos(th).*sin(th)./((1-(ecc^2)*(cos(th)).^2).^1.5);
fun = sqrt(rth.^2 + drdth.^2); 
s_in(1) = 0;

for ith =2:size(th,2)
      s_in(ith) = round(trapz(th(1:ith),fun(1:ith)),rs); 
end

L(iter) = round(s_in(end),rs);  
dx(iter) = dx; 
dx_reg(jreg) = dx;
x = 0:dx(iter):L(iter);  
N = size(x,2); 
c_unif = ones(1,N); 
cL_old = c_unif;
phi_sp(iter) = phi_in;
phi_nuc(iter) = round(phi_sp(iter)+pi/2,rn);

cL_arr = zeros(N,size(tarray,2));
cell_len = zeros(1,size(tarray,2)); 
cell_wid = zeros(1,size(tarray,2));
phi_sp_reg = zeros(1,size(tarray,2)); 
xgrid = x;

cL_arr(:,jreg) = c_unif; 
cell_len(jreg) = cell_a(iter);
cell_wid(jreg) = cell_b(iter);
phi_sp_reg(jreg) = phi_sp(iter);
phi_nuc_reg(jreg) = phi_nuc(iter); 

ctot(iter) = sum(cL_old*dx(iter));
ctot_reg(jreg) = sum(cL_arr(:,jreg)*dx_reg(jreg)); 
jreg = 2; 

%%% for some checks
cscale(iter) = cL_old(50);
len_check(iter) = s_in(end);
%%%%% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%%%%%%% EVOLUTION BEFORE TROUND : SHAPE CHANGES + LGN PATTERNING %%%%%%%%%%
%%%% for tests
% konout = 1;  konin = 1; koffout= 0.02; koffin= 0.02;
tic 
while time(iter) < T_total 
    
    %%%%%%%%%%%%%% inhibition range on periphery %%%%%%%%%%%%%%%%%%%%%%%%%%
    if time(iter) < Tround  
     
    %%%% dna disc is circular: distance of a point on periphery should be
    %%%% less than the radius of this disc for inhibited region 
    dist_PtoC = sqrt((shape(th).*cos(th)).^2+(shape(th).*sin(th)).^2);
    dist_fun = dna_R - dist_PtoC; 
    
    %%%% get off and on rates
    [sunique,iunique,ichk] = unique(s_in);
    offrate_on_s = koffout*ones(size(s_in,2),1); 
    offrate_on_s(dist_fun >=0) = koffin; 
    koff = interp1(sunique,offrate_on_s(iunique),x); 

    onrate_on_s = konout*ones(size(s_in,2),1); 
    onrate_on_s(fun >=0) = konin; 
    kon = interp1(sunique,onrate_on_s(iunique),x); 
        
    else  %%%% from the closest distance of a point on periphery to line segment
    %%%% get the (x,y) coordinates of DNA line segment ends
    N1x = ArmR*cos(phi_nuc(iter)); N1y = ArmR*sin(phi_nuc(iter));
    N2x = ArmR*cos(phi_nuc(iter)+pi); N2y = ArmR*sin(phi_nuc(iter)+pi);

    %%%% Call closest_dist function to calculate smallest distance of each point on
    %%%% periphery to the metaphase plate: RanR = metaphase plate width (7um)
    %%%% + inhibition range (4um)
    dist_PtoN = closest_dist(shape,th,[N1x N1y],[N2x,N2y]);
    dist_fun = RanR - dist_PtoN; 

    %%%% get off and on rates
    [sunique,iunique,ichk] = unique(s_in);
    offrate_on_s = koffout*ones(size(s_in,2),1); 
    offrate_on_s(dist_fun >=0) = koffin; 
    koff = interp1(sunique,offrate_on_s(iunique),x); 

    onrate_on_s = konout*ones(size(s_in,2),1); 
    onrate_on_s(fun >=0) = konin; 
    kon = interp1(sunique,onrate_on_s(iunique),x); 


    end
 
    %lD = sqrt(D./koff); 
    if ran_inhibition ==1   
        %%%% If inhibition on then evolve concentration according to reaction diffusion equation
        cL1=cL_old + rk4_new(cL_old,D,tau,koff,dtc,dx(iter));
        cLhalf=cL_old + rk4_new(cL_old,D,tau,koff,0.5*dtc,dx(iter));
        cL2=cLhalf + rk4_new(cLhalf,D,tau,koff,0.5*dtc,dx(iter));

        %%%%%%%% calculation of error and time step update %%%%%%%%%%%%%%%%%%%%
        erru=norm(cL1-cL2);               
        err_step(iter)=erru/tol_adap; 

        if err_step(iter) ~= 0
             fac=1/(err_step(iter)^0.2);
             dtc=min([dtc*fac 10*dtc]);   %%%% time step update
        end
    else
        %%%% Only geometry and rounding impact on concentration, no Ran
        cL1 = cL_old; 
        dtc = dtreg;
        err_step(iter) = 1;  
        %%% It does not mean anything for constant time, but is required
        %%% for the if condition. Set it to its optimal value 1. 
    end
    
    %%%%%%%%% update only concentration and shape if time(iter) < Tround,and 
    %%%%%% concentration, shape and the spindle angle for time(iter) >
    %%%%%% Tround. Also advance time and iteration number.   
    if(err_step(iter) < 1.5)

    %%%%% update concentration
    cL_new = cL1;
    delt(iter) = dtc;
    
    %%%%% calculate torque on spindle if time > Tround
    if time(iter) < Tround
    phi_sp(iter+1) = phi_sp(1); 
    phi_nuc(iter+1) = phi_nuc(1); 
        
    elseif time(iter) >= Tround
    [res_torque,tdens1,tdens2,MT_ang1,MT_ang2] = Torque_static_pinned_Aug18(cell_a(iter),cell_b(iter),th,s_in...
                                        ,dth,dx(iter),rn,rs,phi_sp(iter),len_sp,Lm,x,L(iter),cL_old);
    %%%%%%%%%%%%%%%%%%%%% update angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phi_sp(iter+1) = round(phi_sp(iter) + mu0*res_torque*delt(iter),rn);
    phi_nuc(iter+1) = round(phi_sp(iter+1) + pi/2,rn);  
   
    end
    
    %%%%% Note that shape parameters are being updated after the 
    %%%%% update of iter, while angle and cL are updated
    %%%%% before it.  So after the following two lines (iter+1)  --> iter. 
    
    time(iter+1) = time(iter) + dtc; 
    iter = iter+1;              
      
    
    %%%%% update shape    
    cell_a(iter) = half_cell_len(time(iter));
    cell_b(iter) = half_cell_wid(time(iter));

    %%%% calculate the new s-grid and new system size
    ecc = sqrt(1-(cell_b(iter)^2)/cell_a(iter)/cell_a(iter));     
    shape = @(th) cell_b(iter)./sqrt(1-(ecc^2)*cos(th).^2);  %%%%  define ellipse  
    rth = shape(th);
    drdth = -cell_b(iter)*(ecc^2)*cos(th).*sin(th)./((1-(ecc^2)*(cos(th)).^2).^1.5);
    fun = sqrt(rth.^2 + drdth.^2); 
    s_in(1) = 0;

    for ith =2:size(th,2)
        s_in(ith) = round(trapz(th(1:ith),fun(1:ith)),rs); 
    end

 
    %%%%% obtain L in every loop
    L(iter) = round(s_in(end),rs);  
    dx(iter) = L(iter)/(N-1); 
    x = 0:dx(iter):L(iter);     %%%% calculate new x, keeping same N
    len_check(iter) = x(end);   %%%% This matchs with the Perimeter(cell_a(iter),cell_b(iter)).
    dim =2;
    
    %%%%%% Change cL_new because of shrinking acc to the local conservation
    %%%%%% rule, i.e.  cL_new(i,iter-1)*dx(iter-1) = cLshrink_new(i)*dx(iter)
    if iter > 1
    cL_new = cL_new*L(iter-1)/L(iter); 
    ctot(iter) = sum(cL_new*dx(iter));
    
    %%%% define cscale which should scale as expected for no-ran case, when
    %%%% increase in c comes only due to rounding/shrinking. 
    cscale(iter) = cL_new(50);
    
    if plt_iter ==1
    figure(1)
    subplot(1,3,1)
    plot(x,cL_new,'linewidth',1);
    xlim([0 100])
    ylim([10 65]) 
    title(['time=',num2str(time(iter))])
    
    subplot(1,3,2)
    Der2=(16*circshift(cL_new,-1,dim)+16*circshift(cL_new,1,dim)...
      -circshift(cL_new,-2,dim)-circshift(cL_new,2,dim)-30*cL_new);
      Der2 = Der2/12/dx(iter)/dx(iter);
      plot(x,Der2,'linewidth',1)% %     hold all
    
    plot(x,Der2)
    xlim([0 100])
    ylim([-30 10]) 
    drawnow
    subplot(1,3,3)
    plot(x,koff)
    xlim([0 100])
    ylim([0 0.2]) 

    drawnow
    end 
%     hold all; 
          
    end
        
    if plt_iter ==1
    tpnt = iter;
    conc = cL_new; 
    tim = time(iter);
    Dynamic_plots_running_sim_with_angle;
    end
  
    %%%%%%%%%%%%%%%%%%%%% Running interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (time(iter-1)<tarray(jreg) && tarray(jreg)<=time(iter))
        Dif = abs(time(iter)-time(iter-1));        
        if Dif >= 1e-3 
            t_irreg = [time(iter-1) time(iter)];             
            ck_int(:,1) = cL_old; 
            ck_int(:,2) = cL_new; 

            [Xo,To] = meshgrid(t_irreg,x); 
            [Xn,Tn] = meshgrid(tarray(jreg),x);
            cL_arr(:,jreg) = interp2(Xo,To,ck_int,Xn,Tn); 

            cell_len(jreg) = round(interp1(t_irreg,[cell_a(iter-1) cell_a(iter)],tarray(jreg)),rs);
            cell_wid(jreg) = round(interp1(t_irreg,[cell_b(iter-1) cell_b(iter)],tarray(jreg)),rs);
            L_reg(jreg) = round(interp1(t_irreg,[L(iter-1) L(iter)],tarray(jreg)),rs);
            dx_reg(jreg) =  round(interp1(t_irreg,[dx(iter-1) dx(iter)],tarray(jreg)),rs);
            ctot_reg(jreg) = sum(cL_arr(:,jreg)*dx_reg(jreg)); 
            
            phi_sp_reg(jreg) = round(interp1(t_irreg,[phi_sp(iter-1) phi_sp(iter)],tarray(jreg)),rs);
            phi_nuc_reg(jreg) = round(phi_sp_reg(jreg)+pi/2,rn); 
            
            xgrid = [xgrid; x];            
            else

            cL_arr(:,jreg) = cL_new;             
            cell_len(jreg) = cell_a(iter-1);
            cell_wid(jreg) = cell_b(iter-1);
            L_reg(jreg) = L(iter-1); 
            phi_sp_reg(jreg) = round(phi_sp(iter-1),rn);
            phi_nuc_reg(jreg) = round(phi_sp_reg(jreg)+pi/2,rn); 
            xgrid = [xgrid; x];
            
            dx_reg(jreg) = dx(iter); 
            ctot_reg(jreg) = sum(cL_arr(:,jreg)*dx_reg(jreg)); 
           
        end
          
        if plt_reg ==1
            tpnt = jreg;
            conc = cL_arr(:,jreg-1); 
            tim =  tarray(jreg-1); 
            
            cl_lo = 0; cl_hi = 30;
            if tim < Tround; Dynamic_plots_running_sim;  drawnow
            else; Dynamic_plots_running_sim_with_angle; drawnow
            end
                        
            if plt_movie ==1
                F = getframe(gcf);
                open(mov); 
                writeVideo(mov,F);
            end
        end
        jreg = jreg+1 ;                                     %%% advance jreg
%         [tarray(jreg) jreg]
    end                                          %%%%% end of interpolation
    cL_old = cL_new;                             %%% replace the solution for next iteration    
    
    end                                               %%%%%% end of error loop 

end
runtime = toc; 

if plt_movie ==1
close(mov)  
end  
  
  
%%%%%%% To save the data   
%%%%% save concentration data 
dlmwrite(['Conc_mu0=',num2str(mu0),'Phi_in=',num2str(rad2deg(phi_initial(phi_ind))),'.csv'],cL_arr)
%%%%% save xgrid 
dlmwrite(['Xgrid_mu0=',num2str(mu0),'Phi_in=',num2str(rad2deg(phi_initial(phi_ind))),'.csv'],xgrid)
%%%%% save cell length, cell width, spindle angle and time array 
Time_mat = [tarray' cell_len'  cell_wid' phi_sp_reg'];
dlmwrite(['Time_quants_mu0=',num2str(mu0),'Phi_in=',num2str(rad2deg(phi_initial(phi_ind))),'.csv'],Time_mat)

% dlmwrite(['Conc_data.csv'],cL_arr); 
% dlmwrite('xgrid.csv',xgrid);
% dlmwrite('Sp_angle.csv',[tarray; phi_sp_reg; phi_nuc_reg])
  
  
end
  
  
  
  
  
  
  
  
  
  
  
  
  
  


