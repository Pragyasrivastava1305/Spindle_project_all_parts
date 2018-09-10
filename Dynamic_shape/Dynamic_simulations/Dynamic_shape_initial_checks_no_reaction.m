%%%%% Main code to specifiy shape dynamics and evolve (a) only spindle
%%%%% dynamics or (b) spindle and LGN dynamics both.  The shape dynamics 
%%%%% is specified as a step change in time at time Tround. In main code
%%%%% the mapping of concentration from initial shape to final shape is
%%%%% encoded. To evolve the system, this code calls the function
%%%%% System_dynamics (to simulate coupled dynamics of LGN and spindle, Dyn_choice=1) and
%%%%% Spindle_dynamics (to simulate only spindle dynamics with no Ran,
%%%%% Dyn_choice=0). Another choice is made between no noise and noisy
%%%%% dynamics of spindle (noi_choice =0, or 1). 
%%%% Parameters set from monopolar model
%%%% rounding decimals for angle and spatial resolution.
rn = 8; rs = 1;
load_sys_params_dyn; 
%%% Assign physical parameters.
ArmR = params.dna_half_length;
RanR = params.inhibition;
Lm = params.astral_Mts; 
len_sp = params.spindle_half_length ; 
D = params.Diffusion; 
mu0 = params.mobility; 
konout = params.kon_far; 
konin = params.kon_near; 
koffout = params.koff_far; 
koffin = params.koff_near; 

mov = VideoWriter(['inhibited_reg_movie'])

%%%%  specify time and space 
T_total = 10; 
Tround = 10; 
Tdelay = 5; 

%%%%% angle of spindle at t=10 mins.
phi_10 = deg2rad(70); 

%%%%% grid and arrays
dtreg = 0.5; 
tarray = 0:dtreg: T_total; 
nsize = size(tarray,2); 
dtc = 0.005; tol_adap = 0.001; 
v1 = [14.84, -0.24];
v2 = [1.62, -0.57]; 

%%%%%%%%%%%%%% initial  set up of system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i =1;  jreg =1; dx = 0.01; 
delt(i) = 0.01;
time(i) = 0; 

cell_a(i) = 0.5*(22.83 + v1(1));
cell_b(i) = 0.5*(22.35 + v2(1));

dna_R = 13; 

shape = @(th) cell_a(i)./sqrt((cell_a(i)/cell_b(i))^2 - ((cell_a(i)/cell_b(i))^2-1)*cos(th).^2);
th = round(linspace(0,2*pi,360*10),rn);
dth = round(th(2)-th(1),3);
rth = shape(th);
den = ((cell_a(i)/cell_b(i))^2- ((cell_a(i)/cell_b(i))^2-1)*cos(th).^2).^1.5;
drdth = -cell_a(i)*((cell_a(i)/cell_b(i))^2 -1)*cos(th).*sin(th)./den;
fun = sqrt(rth.^2 + drdth.^2); 
s_in(1) = 0;

for ith =2:size(th,2)
      s_in(ith) = round(trapz(th(1:ith),fun(1:ith)),rs); 
end

L(i) = round(s_in(end),rs);  
dx(i) = dx; 
x = 0:dx(i):L(i);  
N = size(x,2); 
c_unif = 50*ones(1,N); 
c_array = zeros(N,size(tarray,2));
c_array(:,i) = c_unif;


%%%%%%%%%%%%%%%%%%%%%% BEFORE TROUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Upto the time Tround (time upto which DNA remains round) 
%%%%%%%% Dynamics only consists of shape and patterning, dealt with here. 
% plot(xunit,yunit); axis equal; drawnow
% xlim([-40 40]); ylim([-40 40])

for i =1:size(tarray,2)
    
    
cell_a(i) = 0.5*(22.83 + v1(1)*exp(v1(2)*tarray(i)));
cell_b(i) = 0.5*(22.35 + v2(1)*exp(v2(2)*tarray(i)));

shape = @(th) cell_a(i)./sqrt((cell_a(i)/cell_b(i))^2 - ((cell_a(i)/cell_b(i))^2-1)*cos(th).^2);
th = round(linspace(0,2*pi,360*10),rn);
dth = round(th(2)-th(1),3);
rth = shape(th);
drdth = -cell_a(i)*((cell_a(i)/cell_b(i))^2 -1)*cos(th).*sin(th)./(((cell_a(i)/cell_b(i))^2 - ((cell_a(i)/cell_b(i))^2-1)*cos(th).^2).^1.5);
fun = sqrt(rth.^2 + drdth.^2); 
s_in(1) = 0;

for ith =2:size(th,2)
      s_in(ith) = round(trapz(th(1:ith),fun(1:ith)),rs); 
end

%%% do not change this discretisation. Increasing the resolution needs
%%% rounding place rn to be larger. 
xunit= shape(th).* cos(th) ;
yunit = shape(th).* sin(th) ;

% %%%% plot ellipse as check if needed
% plot(xunit,yunit);  axis equal  
% xlim([-40 40])
% ylim([-40 40])
% drawnow

%%%%% obtain L in every loop
L(i) = round(s_in(end),rs);  
dx(i) = L(i)/(N-1); 
x = 0:dx(i):L(i); 
if i>1
c_array(:,i) = c_array(:,i-1)*L(i-1)/L(i);
end

c_tot(i) = sum(c_array(:,i)*dx(i));
%%%%%  Find inhibition range by finding the interesection of Ran circle
%%%%%  with the ellipse. 

%%%%% calculate the inhibited parts of periphery using phi_nuc(i)
figure(1)
[spr11,spr12,kon] = Rates_dynamic_circular_DNA_May18(x,th,s_in,cell_a(i),cell_b(i),dna_R,konout,konin); 
[spr21,spr22,koff] = Rates_dynamic_circular_DNA_May18(x,th,s_in,cell_a(i),cell_b(i),dna_R,koffout,koffin,'flag'); 

axis equal  
xlim([-40 40])
ylim([-40 40])
title(['time=',num2str(tarray(i))])
grid on
drawnow; 


F = getframe(gcf);
open(mov); 
writeVideo(mov,F);

end
close(mov)  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


