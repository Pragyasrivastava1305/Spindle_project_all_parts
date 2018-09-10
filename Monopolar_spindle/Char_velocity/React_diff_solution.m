function [ctilde,f] =React_diff_solution(dir,num,sp_len,l,D,par1,dx,varargin)
%%%%% This function solves for \tilde{c}(s,t) as defined in the manuscript,
%%%%% which is then used to calulate force on the monopolar spindle. l,D,
%%%%% par1 in this function are to be inputted from the result of
%%%%% optimisation routine. 

%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% dir = the directory containing cell data
%%% num = cell number
%%% sp_len = spindle length for cell num
%%% l = inhibition range
%%% D = diffusion constant
%%% par1 = [koffout koffin konin alpha_i]
%%% varargin is passed as 'flag' if figure needs to be generated. 
  
   
%%%% For troubleshooting. Comment the function line and then 
%%%% input the values of each parameter within the code, to run it as 
% %%%%% a separate program. 
%    dir = 'Cell1'; num=1;     sp_len =13.5; 
%    alpha=50; D=2; l = 10.83; 
%    konout = 1;   konin= 1;
%    koffout = 0.01;  koffin = 0.05;
% %    
% l=11; D=0.1; koffout=0.015; koffin=0.079; 
%   konin=1; konout=1; alpha=49; num=1; sp_len=13.5;
%   dir = 'Cell1';
% % % %%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
lgn_data = load(fullfile(dir,['Cell',num2str(num),'_conc.dat']));
dna_data = load(fullfile(dir,['Cell',num2str(num),'_Mech.dat'])); 
   
xlab = lgn_init(:,1);
% dx=0.5; 
L=floor(xlab(end)/dx)*dx; 
x=(0:dx:L); N=size(x,2); 

Nstep = size(lgn_data,1);                  
for ifr=1:Nstep
    %lgn_aver(ifr,:) = movmean(lgn_data(ifr,:),10);
    lgn_exp(ifr,:) =  interp1(xlab,lgn_data(ifr,:),x);
end

mean_exp = sum(sum(lgn_exp))/L/dna_data(end,1);
lgn_norm = lgn_exp/mean_exp;
   
tvec = dna_data(:,1);    
cen_pos = dna_data(:,2);
angle = dna_data(:,4); 
xNexp=cen_pos+sp_len*cos(pi/180*angle);
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

for iout=1:Nstep   
cL_arr(:,iout) = cL_old;
xN=xNexp(iout);
   
%      xN=36;
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
%   plot(delt,'o-.'); drawnow
% test_t(iout) = sum(delt); 
% clear time; clear delt; 
end
 


% plot(x,cL_new,'-.r','linewidth',2,'MarkerSize',4)
exp_kym = lgn_norm'; 

mean_sim = sum(sum(cL_arr))/L/dna_data(end,1);
sim_kym=cL_arr/mean_sim; 
f=sum(sum((exp_kym-sim_kym).^2));

ctilde = cL_arr; 
T = tvec; 
X = x;
%%%%%%% plot kymograph if flag is passed in varargin
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
