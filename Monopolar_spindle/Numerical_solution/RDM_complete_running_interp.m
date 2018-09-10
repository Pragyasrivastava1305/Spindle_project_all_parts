clear all; close all;
%%%% This code numerically integrates the coupled system of diffusion-reaction
%%%% equation and the force balance on centrosome and nucleus. Initial
%%%% conditions are chosen from the experimental data for one of the cells.
%%%% The code uses adaptive time step combined with RK4. Spatial
%%%% derivatives are approximated by centered finite-differences. The
%%%% limitation on dx is set by the ratio of Diffusion constant and the
%%%% velocity of moving centrosome, which is unknown apriori but is
%%%% primarily set by the inverse of the time scale taua. The result is
%%%% saved only at much fewer times points regularly spaced. Interpolation
%%%% on the go has been used to get solution at regular time points. 
%%%% Spring force is calculated from linear coordinates and uses relative
%%%% position of centrosome and nucleus. Active force is calculated using 
%%%% the periodic coordinate of centrosome. Unlike codes that make
%%%% comparison to experiments,  two loops are not
%%%% needed as position of nucleus is not being updated from the data.
icell = 1;  
plt_chks =0;    %%% set plt_chks to 1 to generate plots while the code runs.
opt_mov = 0; 
dir = (['Cell',num2str(icell)]); 
lgn_init = load(fullfile(dir,['Cell',num2str(icell),'_lgn_init.dat']));
xlab=lgn_init(:,1); 

%%%%%%% system and corrds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx=0.001; L=floor(xlab(end)/dx)*dx;  x=(0:dx:L); N=size(x,2); 
r = 2;                               %%% rounding place, should be set so that 
                                     %%% dx can resolve last decimal place.
                                     %%% dx= 0.1, r=1, dx=0.01, r=2.                                          
dx_nuc=dx;                           %%% dx_nuc <=dx.
tmax=200; dtc=0.01;                  %%% initial choice of dtc 

dtreg=1;  tol=0.001; 
dtm=dtc;                %%%%%% choice for the time step of  mechanical part

%%%%%%%%% Physical parameters of the system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=11;  D=0.01;       alpha = 17.65;                    % inhibition range and diffusion
konout=1; konin=1;  koffout=0.025; koffin=0.129;       % in and out rates
kon=zeros(1,N); koff=zeros(1,N);
tau = 1/koffin; 
%%% for above values of parameters dx = 0.001 works.
dim=2;

Lm=20;                               %%%% fixed length of Microtubule reach
tauc=0;                  %%%% time scale related to rigidity of spring
va=0.9;                  %%%%  Mobility obtained from the mechanical fit = 1.8/2.07
option=1;                           %%% for pdf choose 1 for 1-cdf choose 2

%%%%%%%%%%% Load the experimental data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Exp_data = load(fullfile(dir,['Cell',num2str(icell),'_conc.dat']));
lgn_av = movmean(Exp_data(1,:),10);

cin=interp1(xlab,lgn_av,x);         %%%% initial condition from experiments
cin=alpha*cin/mean(mean(Exp_data));

pos_data = load(fullfile(dir,['Cell',num2str(icell),'_Mech.dat']));
xCin = round(pos_data(1,2),2);           %%% initial position of centrosome
angle = pos_data(1,4); 
xNin=round(xCin+13.5*cos(pi/180*angle),2);  %%% initial position of nucleus
epsilon=abs(mean(13.5*cos(pi/180*pos_data(:,4))));  %%% separation between nucleus and centrosome in 1-d model.
% eps = 13.5;
MT_ends = load('Prob_dist_MTends.dat');
MT_x=MT_ends(:,1);    MT_pdf=MT_ends(:,2); 

% mov=VideoWriter(['D=',num2str(D),'va=',num2str(10*va),'.avi']);

%%%%%% Regular time array 
t_reg = 0:dtreg:tmax; 
i=1;  time(i)=0;   jreg=2;
delt(i) = dtc;
%%%% variables used for evolving the code. 
xN_per_old = xNin;   xC_per_old = xCin; 
xN_lin_old = xN_per_old;   xC_lin_old = xC_per_old; 
cL_old=cin;

%%%%Variables used for storing data at regular time points. Initialiase
%%%% for jreg=1.
xN_per_reg(1)=xNin;   xC_per_reg(1)=xCin; 
xN_lin_reg(1)=xNin;   xC_lin_reg(1)=xCin;
cL_arr=zeros(N,size(t_reg,2));
cL_arr(:,1)=cL_old;

mean_pos(1) = 0.5*(xC_lin_reg(1)+xN_lin_reg(1)); 
rel_pos(1) = xC_lin_reg(1) - xN_lin_reg(1); 

vmean(1) = 0;
vrel(1) = 0; 
tic

while (time(i)<tmax)   
    
    %%%%% advance nucleus and centrosome using concentration cL_old
    if (xC_lin_old>=xN_lin_old)
        fspr(i)=(xC_lin_old-xN_lin_old-epsilon);
    elseif (xC_lin_old<xN_lin_old)
        fspr(i)=(xC_lin_old-xN_lin_old+epsilon);
    end
  
    [fact(i),Ipos,Ineg,xper,Pdfper,cse] = Fact_periodic_Mtdist(...
                                 xC_per_old,MT_x,MT_pdf,dx_nuc,x,dx,r,cL_old);
        
    %%%%% update linear position
    if tauc ~=0
    xN_lin_new=xN_lin_old+dtc*fspr(i)/tauc; 
    xC_lin_new=xC_lin_old+dtc*(va*fact(i)-fspr(i)/tauc); 
    
    elseif tauc ==0    %%%%  rigid spring case
    xC_lin_new = xC_lin_old + dtc*va*fact(i) ;  
    xN_lin_new = xC_lin_old -epsilon; 
    end

    %%%% for checks
    mean_pos(i+1) = 0.5*(xC_lin_new+xN_lin_new); 
    rel_pos(i+1) = xC_lin_new - xN_lin_new; 
    vmean(i+1) = (mean_pos(i+1) - mean_pos(i))/dtc;
    vrel(i+1) = (rel_pos(i+1) -rel_pos(i))/dtc;

        if (mod(i,500)==0 && plt_chks ==1) 
%          cse
         time(i);
         figure(1)
         %%% mean velocity and pulling forces   
         subplot(1,3,1)
         hold off        
         if tauc ~=0 
         plot(time(1:end),va*fact/2,'bo','linewidth',1.5,'MarkerSize',2); hold on
         else 
         plot(time(1:end),va*fact,'bo','linewidth',1.5,'MarkerSize',2); hold on
         end
         plot(time(1:end),vmean(2:end),'k','linewidth',0.8,'MarkerSize',2); 
         xlim([0 tmax]);  
         axis square; 
         xlabel(['t(min)']);  
         title([' time=',num2str(time(i))])
         ylabel(['Velocity of mean position'])
         drawnow  
         
         %%% mean and relative velocities   
         subplot(1,3,2)
         hold off
         if tauc ~=0
         plot(time(1:end),(va*fact)-2*(fspr/tauc),'ro','linewidth',1.5,'MarkerSize',2); hold on
         else
         plot(time(1:end),zeros(1,size(time(1:end),2)),'ro','linewidth',1.5,'MarkerSize',2); hold on 
         end
         
         plot(time(1:end),vrel(2:end),'k','linewidth',0.8,'MarkerSize',2); 
         xlim([0 tmax]);  
         axis square; 

         xlabel(['t(min)']); 
         ylabel(['Relative velocity'])
         title([' time=',num2str(time(i))])
         drawnow 
         
         subplot(1,3,3)
         plot(x,cL_old,'linewidth',2)
         axis square;
%             subplot(1,3,3)
%             plot(time(1:end-1),delt)
%             xlabel(['t(min)']); 
%             ylabel(['\Delta t'])
%             axis square;    
        end
   
    %%%%%%% update periodic coordinates from linear coordinates.
    fac_N=xN_lin_new/L;   
    %%%% integer part of fac_N is the number of periodic image. 
 
    if fac_N>=0
        int_N=floor(fac_N); del_N=fac_N-int_N;
        xN_per_new=del_N*L;
        elseif fac_N<0
        int_N=ceil(fac_N);del_N=fac_N-int_N;
        xN_per_new=L*(1+del_N);
    end

    fac_C=xC_lin_new/L;
 
    if fac_C>=0
        int_C=floor(fac_C); del_C=fac_C-int_C;
        xC_per_new=del_C*L;
    elseif fac_C<0
        int_C=ceil(fac_C); del_C=fac_C-int_C;
        xC_per_new=L*(1+del_C);
    end
  
   
    %%%%%% update kon according to the nucleus position
    kon=Rates(konin,konout,x,dx,xN_per_old,l,L); 
    koff=Rates(koffin,koffout,x,dx,xN_per_old,l,L);
    %lD = sqrt(D./koff);
 
    %%%%%% chemical part : calculated from periodic coordinates of DNA
    %%%%%% error calculation for adaptive step:rk4 in time,explicit in space
    cL1=cL_old+rk4_new(cL_old,D,tau,koff,dtc,dx);
    cLhalf=cL_old+rk4_new(cL_old,D,tau,koff,0.5*dtc,dx);
    cL2=cLhalf+rk4_new(cLhalf,D,tau,koff,0.5*dtc,dx);
   
    erru=norm(cL1-cL2);               
    error(i)=erru/tol;    
    if error(i)~=0
       fac=1/(error(i)^0.2);
       dtc=min([dtc*fac 10*dtc]);   %%%% time step update
    end
   
    %%%%% iterate solution if error is sufficiently small %%%%%%%%%%%%%%
    %%%%% and perform interpolation
    if(error(i)<1.5)
%         cL_new=cL_old+rk4_new(cL_old,D,tau,koff,dtc,dx);
        cL_new = cL1; 
        delt(i)=dtc;       
        time(i+1)=time(i)+dtc;
    
        %%%% check if t_reg(jreg) falls between time(i) and time(i-1)
        if (time(i)<t_reg(jreg) && t_reg(jreg)<=time(i+1))
        t_reg(jreg);
        Dif = abs(time(i+1)-time(i)); 
            
        %%%% if (time(i+1)-time(i)) is too small, then instead of
        %%%% interpolating, we allocate the current step in regular
        %%%% array by the values of variables at time(i+1);             
        if Dif >= 1e-3 
        t_irreg = [time(i) time(i+1)]; 

        ck_int(:,1) = cL_old; 
        ck_int(:,2) = cL_new; 

        [Xo,To] = meshgrid(t_irreg,x);
        [Xn,Tn] = meshgrid(t_reg(jreg),x); 
        cL_vec = interp2(Xo,To,ck_int,Xn,Tn);
        cL_arr(:,jreg) =  cL_vec;
                           
        xN_lin_reg(jreg) = interp1(t_irreg,[xN_lin_old xN_lin_new],t_reg(jreg)); 
        xC_lin_reg(jreg) = interp1(t_irreg,[xC_lin_old xC_lin_new],t_reg(jreg)); 
        xN_per_reg(jreg) = interp1(t_irreg,[xN_per_old xN_per_new],t_reg(jreg)); 
        xC_per_reg(jreg) = interp1(t_irreg,[xC_per_old xC_per_new],t_reg(jreg)); 

        clear ck_int
        else
        cL_arr(:,jreg) = cL_new;
        xN_lin_reg(jreg) = xN_lin_new; 
        xC_lin_reg(jreg) = xC_lin_new; 
        xN_per_reg(jreg) = xN_per_new; 
        xC_per_reg(jreg) = xC_per_new; 
        end
            
        %%%% forces at regular time points
        [fact_reg(jreg),Ipos,Ineg,xper,Pdfper,cse] = Fact_periodic_Mtdist(...
                                 xC_per_reg(jreg),MT_x,MT_pdf,dx_nuc,x,dx,r,cL_old);
            
        if (xC_lin_reg(jreg)>=xN_lin_reg(jreg))
        fspr_reg(jreg)=(xC_lin_reg(jreg)-xN_lin_reg(jreg)-epsilon);
        elseif (xC_lin_reg(jreg)<xN_lin_reg(jreg))
        fspr_reg(jreg)=(xC_lin_reg(jreg)-xN_lin_reg(jreg)+epsilon);
        end
            
        %%% advance jreg
           
%         %%%%%% make plots 
%         figure(1)
%             
%         Plots_MTdist_running_interp; 
%         %%% capture frame in movie
%         
%         F = getframe(gcf);
%         open(mov);
%         writeVideo(mov,F);
            
        jreg = jreg+1 
%         plot(cL_new); drawnow

        end
        
        %%% advance i
        i=i+1;
        %%%%%% replace old with new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cL_old=cL_new;
        xN_lin_old = xN_lin_new; 
        xC_lin_old = xC_lin_new; 
        xN_per_old = xN_per_new; 
        xC_per_old = xC_per_new; 
          
   
    end  
   
    runtime = toc;
end
close(fig)
close(mov);



 
%  Mat_N=zeros(size(xN_lin_reg,2),2);
%  Mat_C=zeros(size(xC_lin_reg,2),2);
%  Mat_N(:,1)=xN_lin_reg;  Mat_N(:,2)=xN_per_reg;
%  Mat_C(:,1)=xC_lin_reg;  Mat_C(:,2)=xC_per_reg;
%  
%  dlmwrite(fullfile(dir,['cL_dx=',num2str(100*dx),'.dat']),cL_arr,'delimiter','\t')
%  dlmwrite(fullfile(dir,['f_act=',num2str(100*dx),'.dat']),fact_reg,'delimiter','\t')
%  dlmwrite(fullfile(dir,['f_spr=',num2str(100*dx),'.dat']),fspr_reg,'delimiter','\t')
%  dlmwrite(fullfile(dir,['xN_pos=',num2str(100*dx),'.dat']),Mat_N,'delimiter','\t')
%  dlmwrite(fullfile(dir,['xC_pos=',num2str(100*dx),'.dat']),Mat_C,'delimiter','\t')
 
save(['D=',num2str(D),'_va=',num2str(va),'_tauc=',num2str(tauc),'.mat'])
 
 
 

 
 
 
 
 
 
 
 
 
 






