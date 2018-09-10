%%% Evolution of LGN concentration and spindle angle using rk4-adpative, with a 
%%% time delay introduced between inhibition and spindle rotation
clear all; 
%%%% rounding decimals for angle and spatial resolution.
rn = 8; rs = 1;

%%% ASSIGN PHYSICAL PARAMETERS OF THE SYSTEM
ArmR = 9;    %%% half length of metaphase plate
RanR = 11;   %%% Range of RAN inhibiton
Lm = 15;     %%% Length of MTs
len_sp = 7;  %%% half length of spindle
D = 0.01;    %%% Diffusion constant
mu0 = 0.00063; %%% mobility
Taud = 9;     %%% Delay between NEP and force application on spindle
konout = 1; konin = 1;          %%% On rates outside and inside inhibiton range
koffout = 0.039; koffin = 0.172;  %%% off rates outside and inside. 

%%%%%% vector of initial orientation of spindle
phi_initial = round(linspace(0,pi/2,45),rn); 
Asp_vec = 1.2:0.2:5; 
wid_vec = 12.5; 

%%%%% major and minor axes of ellipse
a = 20; b=10;

%%%% system specification
%%% Time
tmax = 100;
dtc = 0.01; 
dtreg = 1; t_reg = 0:dtreg:tmax;
nsize = size(t_reg,2);
tol_adap = 0.001; 

%%% Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx= 0.05;
%%%%% specify shape of ellipse, calculate perimeter=linear system size.  
ecc = (1-b^2/a/a); 
shape = @(th) b./sqrt(1 - (ecc^2)*cos(th).^2);
th = round(linspace(0,2*pi,360*10),rn);
dth = round(th(2)-th(1),3);

rth = shape(th);
drdth = -a*((a/b)^2 -1)*cos(th).*sin(th)./(((a/b)^2 - ((a/b)^2-1)*cos(th).^2).^1.5);
fun = sqrt(rth.^2 + drdth.^2); 
s(1) = 0;

for ith =2:size(th,2)
      s(ith) = round(trapz(th(1:ith),fun(1:ith)),rs); 
end
[sunique,iunique,jch] = unique(s);  %%%% unique values in s needed for a later interpolation 
%%% do not change this discretisation. Increasing the resolution needs
%%% rounding place rn to be larger. 
xunit = shape(th).* cos(th) ;
yunit = shape(th).* sin(th) ;

%%%% plot ellipse as check if needed
% plot(xunit,yunit); axis equal

%%%%%% system and coords %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% rounding place of L and dx and s should be same. s rounding is enters
%%%% the function to calculate torque as rs. x vector has different length
%%%% than s. s has non-uniform grid while x has uniform grid. 
L = round(s(end),rs);  
x = 0:dx:L; 
N = size(x,2); 

%%%% define the 0 of torque for given a and b. Its a small discretisation
%%%% error. 
c_unif = ones(1,size(x,2)); 
tol_torque =  abs(Torque_static_pinned_Aug18(a,b,th,s,dth,dx,rn,rs,0,len_sp,Lm,x,L,c_unif));
iphi = 10; 
%for iphi = 1:size(phi_initial,2)
% iphi=45;
%[iWid iAsp iphi]

%%%%% initial conditions : i keeps track of adaptive index, jreg is
%%%%% index for regular time grid. 
i = 1;     jreg =1;  
time(i) = 0;
cL_old =  (ones(1,N)); 

phi_sp(i) =  round(phi_initial(iphi),rn); 
phi_nuc(i) = round(phi_sp(i)+pi/2,rn);    

cL_arr = zeros(N,nsize) ;
phi_sp_reg = zeros(1,nsize); 

phi_sp_reg(jreg) = phi_sp(1); 
cL_arr(:,jreg) = cL_old;
jreg = jreg+1;    %%%% jreg index is +1'd because it is the second point that needs the interpolation.
                  %%%% @ jreg = 1 cL_arr = cL_old; 
tic

while time(i) < tmax

%%%% get the (x,y) coordinates of DNA line segment ends
N1x = ArmR*cos(phi_nuc(i)); N1y = ArmR*sin(phi_nuc(i));
N2x = ArmR*cos(phi_nuc(i)+pi); N2y = ArmR*sin(phi_nuc(i)+pi);

%%%% Call closest_dist function to calculate smallest distance of each point on
%%%% periphery to the metaphase plate
dist_PtoN = closest_dist(shape,th,[N1x N1y],[N2x,N2y]);
fun = RanR - dist_PtoN; 

%%%% get off and on rates
offrate_on_s = koffout*ones(size(s,2),1); 
offrate_on_s(fun >=0) = koffin; 
koff = interp1(sunique,offrate_on_s(iunique),x); 

onrate_on_s = konout*ones(size(s,2),1); 
onrate_on_s(fun >=0) = konin; 
kon = interp1(sunique,onrate_on_s(iunique),x); 

lD = sqrt(D./koff); 
cL1 = cL_old + rk4_RDM(cL_old,lD,kon,koff,dtc,dx);
cLhalf = cL_old + rk4_RDM(cL_old,lD,kon,koff,0.5*dtc,dx);
cL2 = cLhalf + rk4_RDM(cLhalf,lD,kon,koff,0.5*dtc,dx);

erru=norm(cL1-cL2);               
err_step(i)=erru/tol_adap; 

if err_step(i) ~= 0
    fac=1/(err_step(i)^0.2);
    dtc=min([dtc*fac 10*dtc]);   %%%% time step update
end
    
%%% update the solution for sufficiently small error
if(err_step(i) < 1.5)
cL_new = cL1;
delt(i) = dtc;
time(i+1) = time(i) + dtc; 
 
[res_torque,tdens1,tdens2,MT_ang1,MT_ang2] = Torque_static_pinned_Aug18(a,b,th,s...
                                        ,dth,dx,rn,rs,phi_sp(i),len_sp,Lm,x,L,cL_old);
                         
%%% set res_torque to 0 if time < Taud
if time(i) < Taud
 res_torque = 0;  
end

%%% update angle
%%% phi_sp(i+1) = round(phi_sp(i)+5*delt(i),rn); 
phi_sp(i+1) = round(phi_sp(i) + mu0*res_torque*delt(i),rn);
phi_nuc(i+1) = round(phi_sp(i+1) + pi/2,rn); 

%%%%%%%%%%%%%%%%%%%%%%% Running interpolation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% if time(i+1) - time(i) is too small, then instead of interpolating, allocate
%%%% the current step in regular array by the values of variables at time(i+1).

%%%% check if t_reg(jreg) falls between time(i) and time(i+1)
if (time(i)<t_reg(jreg) && t_reg(jreg)<=time(i+1))
    % [t_reg(jreg) time(i+1)]
    Dif = abs(time(i+1)-time(i)); 
           
    if Dif >= 1e-3 
        
    t_irreg = [time(i) time(i+1)]; 
    ck_int(:,1) = cL_old; 
    ck_int(:,2) = cL_new; 
    [Xo,To] = meshgrid(t_irreg,x);
    [Xn,Tn] = meshgrid(t_reg(jreg),x); 

    %%% use interp1 to interpolate on time at each space grid.
    cL_arr(:,jreg) = interp2(Xo,To,ck_int,Xn,Tn);
    phi_sp_reg(jreg) = round(interp1(t_irreg,[phi_sp(i) phi_sp(i+1)],t_reg(jreg)),rn);
    phi_nuc_reg(jreg) = mod(phi_sp_reg(jreg) + pi/2,2*pi); 
  
   % plot(x,cL_new); drawnow
    else
        
    cL_arr(:,jreg) = cL_new; 
    phi_sp_reg(jreg) = phi_sp(i+1);  
    phi_nuc_reg(jreg) = mod(phi_sp_reg(jreg) + pi/2,2*pi); 
  
    %plot(x,cL_new); drawnow
    end                                          %%%% end of time diff check condition.

    
iframe = jreg;   
plot(x,cL_new); drawnow
%Plots_movie_Static; drawnow
jreg = jreg+1;                                   %%% advance jreg

end                                              %%%%% end of interpolation

%%% time(i);   
i=i+1;                                            %%%% advance i
cL_old = cL_new;                                  %%%% replace the solution

end                                               %%%%%% end of error loop 

    
end
runtime =  toc; %%%% end of time loop.


%%%% to save data : use gridspacing larger than used for simulation. 
xsave = 0:0.5:L;                    %%%% save data at 0.1 um distance
[X0,T0] = meshgrid(x,t_reg); 
[XS,TS] = meshgrid(xsave,t_reg); 
cL_save = interp2(X0,T0,cL_arr',XS,TS);

NS = size(xsave,2); 
save_data = zeros(NS+1,nsize);
save_data(1:NS,:) = cL_save';
save_data(end,:)  = phi_sp_reg;

% %%% save data in a different directory specified by dir
% dlmwrite(fullfile(dir,['iPhi=',num2str(iphi),'RanR=',num2str(RanR),'_Taud=',num2str(Taud)...
%     ,'_=a',num2str(a),'_=b',num2str(b),'_Lm=',num2str(Lm),'.csv']),save_data,'delimiter','\t'); 
%  
% save(fullfile(dir,['iPhi=',num2str(iphi),'RanR=',num2str(RanR),'_Taud=',num2str(Taud)...
%     ,'_=a',num2str(a),'_=b',num2str(b),'_Lm=',num2str(Lm),'.mat'])); 


% %%%%%%% save data in same directory
%  dlmwrite(['iPhi=',num2str(iphi),'RanR=',num2str(RanR),'_Taud=',num2str(Taud)...
%      ,'_A=',num2str(Asp),'_b=',num2str(b),'_Lm=',num2str(Lm),'.csv'],save_data,'delimiter','\t'); 
%  
% %  save(['iPhi=',num2str(iphi),'RanR=',num2str(RanR),'_Taud=',num2str(Taud)...
% %      ,'_=a',num2str(a),'_=b',num2str(b),'_Lm=',num2str(Lm),'.mat']); 
% clear time; 

%end

%clear x ck_int cL_save;
% end
% end





