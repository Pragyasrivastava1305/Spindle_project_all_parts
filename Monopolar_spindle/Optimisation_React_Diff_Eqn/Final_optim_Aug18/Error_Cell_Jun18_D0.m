function f = Error_Cell_Jun18_D0(dir,num,sp_len,l,par1,varargin)
%%%% This code is same as Error_Cell_Jan18_single, except that it does not
%%%% use implicit scheme for solution updates in time. This code deals only
%%%% with D = 0 case where the reaction diffusion terms has only local
%%%% terms. The implicit method, if employed here, will make delt increase
%%%% by its maximum allowed value, since error between solution updated by
%%%% single time-step and the solution updated by two half time-steps is 0.

%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 konout = 1; koffout = par1(1);  
 koffin = par1(2);  konin=par1(3); alpha=par1(4);
%  l=11; D=0; koffout=0.015; koffin=0.079; 
%  konin=1; konout=1; alpha=49; num=1; sp_len=13.5;
%  dir = 'Cell1';

dtreg = 3;
clear lgn_aver lgn_exp 
lgn_init = load(fullfile(dir,['Cell',num2str(num),'_lgn_init.dat'])); 
lgn_data = load(fullfile(dir,['Cell',num2str(num),'_conc.dat']));
dna_data = load(fullfile(dir,['Cell',num2str(num),'_Mech.dat'])); 
   
xlab = lgn_init(:,1);
dx=0.1; 
L=floor(xlab(end)/dx)*dx; 
x=(0:dx:L); N=size(x,2); 

Nstep = size(lgn_data,1);                  
for ifr=1:Nstep
    lgn_aver(ifr,:) = movmean(lgn_data(ifr,:),10);
    lgn_exp(ifr,:) =  interp1(xlab,lgn_aver(ifr,:),x);
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

%%%%% initial condition of the simulation 

cin = alpha*lgn_norm(1,:); 
cL_arr = zeros(N,Nstep); 
cL_arr(:,1) = cin; 

for iout = 1:Nstep-1
    %%%%%% update kon according to the nucleus position
    kon = Rates(konin,konout,x,dx,xNexp(iout),l,L); 
    koff = Rates(koffin,koffout,x,dx,xNexp(iout),l,L);
 
    %%%%%% chemical part : calculated from periodic coordinates of DNA
    %%%%%% error calculation for adaptive step:rk4 in time,explicit in space
 
    RHS(:,iout) = kon' - koff'.*cL_arr(:,iout); 
    cL_arr(:,iout+1)=cL_arr(:,iout)+dtreg*RHS(:,iout);
   
end

%%%%  to make cL_arr non-dimensional in the same way as in supplementary. 
cL_arr = koffin*cL_arr; 
exp_kym = lgn_norm'; 

mean_sim = sum(sum(cL_arr))/L/dna_data(end,1);
sim_kym=cL_arr/mean_sim; 
f=sum(sum((exp_kym-sim_kym).^2));


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
 
 
 
 
 
 
 






