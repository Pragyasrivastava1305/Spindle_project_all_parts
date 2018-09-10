function f=Error_Cell_Jun18_single(dir,num,sp_len,l,D,par1,dx,varargin)

%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% dir = the directory containing cell data
%%% num = cell number
%%% sp_len = spindle length for cell num
%%% l = inhibition range
%%% D = diffusion constant
%%% par1 = [koffout, koffin, konin, alpha_num]
%%% varargin is passed as 'flag' if figure needs to be generated. 
  
%%%% For troubleshooting. Comment the function line and then 
%%%% input the values of each parameter within the code, to run it as 
%%%% a separate program. 

%   num=1; sp_len=13.5;
%   dir = 'Cell1';
%   l=11; D=1;
%   dx =0.2;
% 
% % %%%%% most recent optimisation results : Cell 1
%  koffout=0.07; koffin=0.3; alpha=5;
%  konin=1; konout=1;

% %%%%% Parameters from GS
%  koffout = 0.033;  koffin = 0.2; alpha= 87;
%  koin = 1; konout =1;
  
%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
konout = 1; koffout = par1(1);  
koffin = par1(2);  konin=par1(3); alpha=par1(4);
   
tmax = 3;            % time between two frames. 
dtc = 0.01;          % initial time step.
tol = 0.001;         % tolerance for adaptive step size
 
%%%%%%%%%%% IMPORT EXPERIMENTAL DATA FOR num %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  a) Load initial conditions :  define system coords from first column
%%%  and define the numerical grid. b) Load complete kymograph :  
%%%  interpolate concentration at numerical grid points c) Load data for
%%% positions of centrosome and spindle angle. Define and periodise 
%%% nucleus position
 
clear lgn_aver lgn_exp 
lgn_init = load(fullfile(dir,['Cell',num2str(num),'_lgn_init.dat'])); 

%%%%% lgn_init = [experimental space grid, initial lgn concentration]
%%%%% space grid contains grid points around periphery in experiments.
%%%%% Initial lgn concentration = Initial lgn on the above grid
%%%%% points hence in units of intensity.

lgn_data = load(fullfile(dir,['Cell',num2str(num),'_conc.dat']));
%%%% lgn_data = concentration map lgn_data(i,j), i runs over experimental
%%%% time points, j over grid point as recorded in lgn_init(:,1). concentration is 
%%%% in units of intensity

dna_data = load(fullfile(dir,['Cell',num2str(num),'_Mech.dat'])); 
%%%% dna_data =[time, periodic_centrosome_position, unfolded_centrosome_displacement, monopolar angle]
%%%% DNA position is calculated from periodic centrosome_position and
%%%% monopolar angle   
xlab = lgn_init(:,1);         %%%  experimental grid
L=floor(xlab(end)/dx)*dx; 

%%%  system size. In experiments the positions are measured upto 4th
%%%  decimal point by ImageJ. In the code we round up based on the
%%%  numerical grid size. An example is : if xlab(end)=11.5670 and dx=0.5
%%% xlab(end)/dx = 23.140, floor(xlab(end)/dx)= 23. Always the lower
%%% integer is chosen. Chosing the higher itneger (e.g. 24 in this example will 
%%% lead to position which are bigger than the periphery). 
%%% Finally floor(xlab(end)/dx)*dx= 23*dx = 11.5. Thus, this way of defining
%%% results in 11.5670 ~ 11.5.

x=(0:dx:L);         %%%%% numerical grid
N=size(x,2);     

Nstep = size(lgn_data,1);                  
for ifr=1:Nstep
    %%% interpolate experimental LGN concentration over numerical grid.
    lgn_exp(ifr,:) =  interp1(xlab,lgn_data(ifr,:),x);  
end

%%%% experimental spatiotemporal mean. 
mean_exp = dx*3*sum(sum(lgn_exp))/L/dna_data(end,1);

%%%% experimental concentration normalised by spatiotemporal mean. Dimensionless quantity
lgn_norm = lgn_exp/mean_exp;    

tvec = dna_data(:,1);    
cen_pos = dna_data(:,2);
angle = dna_data(:,4); 

%%%%  calculate DNA position
xNexp=cen_pos+sp_len*cos(pi/180*angle);

%%%% Make DNA position periodic
for ifr=1:Nstep
    if (xNexp(ifr) < xlab(1))
           xNexp(ifr) = xlab(end) - abs(xNexp(ifr));
    elseif (xNexp(ifr) > xlab(end))
           xNexp(ifr) = xNexp(ifr) - xlab(end);
    end
end

%%% Initial condition of simulation. To compare with old manual minimisation 
%%% use cin = lgn_exp(1,:)/11000;
cin = alpha*lgn_norm(1,:); 
   
%%% Solving for LGN concentration in time and space.  
%%% Nucleus position is taken from experimental data and is updated
%%% every 3 minutes. The inner loop between two updates of nucleus
%%% position evolves reaction diffusion for LGN combining 
%%% rk4 and adaptive time step  and euler discretisation in space
cL_old = cin;
cL_arr = zeros(N,Nstep);
tsim(1) = 0;

for iout=1:Nstep   
cL_arr(:,iout) = cL_old;
xN=xNexp(iout);
   
%%%%%% update kon according to the nucleus position
kon = Rates(konin,konout,x,dx,xN,l,L);
koff = Rates(koffin,koffout,x,dx,xN,l,L);
lD=sqrt(D./koff);
i=1;                          %%%%% initialise inner loop for 3 minutes
time(i)=0;   
        
while (time(i) < tmax)                           
cL1=cL_old + rk4_RDM(cL_old,lD,kon,koff,dtc,dx);
            
cLhalf=cL_old + rk4_RDM(cL_old,lD,kon,koff,0.5*dtc,dx);
cL2=cLhalf + rk4_RDM(cLhalf,lD,kon,koff,0.5*dtc,dx);

erru=norm(cL1-cL2);               
error(i)=erru/tol; 

if error(i) ~= 0
    fac=1/(error(i)^0.2);
    dtc=min([dtc*fac 10*dtc]);   %%%% time step update
end
            
if(error(i) < 1.5)
    cL_new = cL1;
    delt(i) = dtc;
    time(i+1) = time(i) + dtc;
    i=i+1;
    cL_old = cL_new;
end

end
 
%%%%% create tsim vector by adding time(end)(~3) to previous tsim point.  
if iout >1
tsim(iout) = tsim(iout-1)+time(end);  
end

%%%%% uncomment to plot cL profile in x for each time
%   subplot(1,2,1)
%   hold off
%   plot(x,cL_new,'-.r','linewidth',2,'MarkerSize',4)
%   hold on 
%   scatter(xN,0); 
%   xlim([0 L])
%   axis square
%   drawnow
%              
%   subplot(1,2,2) 
%   hold off
%   plot(x,kon,'.','linewidth',1.5); hold on; plot(x,koff,'.','linewidth',1.5);
%   scatter(xN,0); 
%   xlim([0 L]); ylim([0 1])
%   axis square
%   drawnow
%   figure; 

%%%%%% uncomment to see the convergence of delt. 
%plot(delt,'o-.'); drawnow

%%%% clear time and delt for the next 3 min interval 
clear time; clear delt; 
end
 
%%%%% interpolate cL_arr(x,tsim) on cL_arr(x,tvec), where tvec is time
%%%%% vector from experiments, taken exactly at 3 mins interval.
[Xo,To] = meshgrid(tsim,x);
[Xn,Tn] = meshgrid(tvec,x); 
cL_arr_reg = interp2(Xo,To,cL_arr,Xn,Tn);  

exp_kym = lgn_norm'; 

mean_sim = dx*3*sum(sum(cL_arr_reg))/L/tvec(end);
sim_kym = cL_arr_reg/mean_sim; 

f = sum(sum((exp_kym-sim_kym).^2));

% %%%%% plot kymograph if flag is passed in varargin
flag = ismember('flag',varargin); 
if flag == true
   figure; 
   set(0,'DefaultAxesFontSize',14,'DefaultTextFontSize',14);
   subplot(2,1,1)
   imagesc(tvec,x,exp_kym); 
   axis xy; axis equal tight; 
   colorbar; caxis([0 max(max(exp_kym))])
   hold on; ylim([0 x(end)])
   plot(tvec,xNexp,'.w','linewidth',2,'MarkerSize',16)
   xlabel('time(min)'); ylabel('x(\mum)')
      
   subplot(2,1,2)    
   imagesc(tvec,x,sim_kym); 
   axis xy; axis equal tight; hold on; 
   colorbar; caxis([0 max(max(exp_kym))])
   plot(tvec,xNexp,'.w','linewidth',2,'MarkerSize',16)
   xlabel('time(min)'); ylabel('x(\mum)')
   drawnow
end
end
