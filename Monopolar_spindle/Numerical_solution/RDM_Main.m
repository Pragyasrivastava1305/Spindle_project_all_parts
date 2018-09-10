clear all; close all;

%%%%%%% system and corrds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=72; dx=0.2; 
x=(0:dx:L);
N=size(x,2); 
dx_nuc=0.01;
tmax=213; Nstep=72; dtc=0.01;                      %%%%initial choice of dtc 
dtreg=1;  tol=0.001; 
dtm=dtc;                %%%%%% choice for the time step of  mechanical part

%%%%%%%%% Physical parameters of the system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=15; lD=sqrt(10);     %%%%%half witdh of inhibited region:nucleus + RanGTP
konin=1; konout=10;  dim=2;
koffout=1; koffin=5;

kon=zeros(1,N); koff=zeros(1,N);
Lm=15;                               %%%% fixed length of Microtubule reach
eps=13.5;                    %%%% separation between centrosome and nucleus
tauc=1;                       %%%% time scale related to rigidity of spring
taua=1.5;
%%%%%%%%%%% Loading the experimental data  %%%%%%%%%%
lgn_init=load('LGN_init_lab.dat'); 
xlab=lgn_init(:,1); lgn_in=lgn_init(:,2); lgn_av=movmean(lgn_in,10);
cin=interp1(xlab,lgn_av,x);        %%%% initial conditions from experiments
  cin=cin/lgn_av(1);

%  cin=rand(1,N);
% xCin=14.5;  xNin=0; 
xCin=30.65; xNin=18.49;

cL_old=cin;
cL_arr=zeros(N,Nstep);

mov=VideoWriter(['RDM_main_exp_init.avi']);

i=1;                          
time(i)=0;   
xN_per(i)=xNin;   xC_per(i)=xCin; 
xN_lin(i)=xN_per(i); xC_lin(i)=xC_per(i);

cL_arr(:,i)=cL_old;

%%%%% calculate initial forces.
if (xC_lin(i)>=xN_lin(i))
     fspr(i)=(xC_lin(i)-xN_lin(i)-eps);
%      [fact(i),Ip,In]=fcen_per(xC_per(i),Lm,x,dx,cL_old); 
     [fact(i),Ip,In]=Fact_periodic(xC_per(i),Lm,x,dx,dx_nuc,cL_old);
elseif (xC_lin(i)<xN_lin(i))
     fspr=(xC_lin(i)-xN_lin(i)+eps);
     [fact(i),Ip,In]=Fact_periodic(xC_per(i),Lm,x,dx,dx_nuc,cL_old);
end
  


while (time(i)<tmax)   
    
 
    
 %%%%%% update kon according to the nucleus position
 kon=Rates(konin,konout,x,dx,xN_per(i),l,L); 
 koff=Rates(koffin,koffout,x,dx,xN_per(i),l,L);
 %%%%%% chemical part : calculated from periodic coordinates of DNA
 %%%%%% error calculation for adaptive step:rk4 in time,explicit in space
 cL1=cL_old+rk4_RDM(cL_old,lD,kon,koff,dtc);
   
 cLhalf=cL_old+rk4_RDM(cL_old,lD,kon,koff,0.5*dtc);
 cL2=cLhalf+rk4_RDM(cLhalf,lD,kon,koff,0.5*dtc);
   
 erru=norm(cL1-cL2);               
 error(i)=erru/tol; 
   
 if error(i)~=0
       fac=1/(error(i)^0.2);
       dtc=min([dtc*fac 10*dtc]);   %%%% time step update
 end
   
 %%%%%%%% iterate solution if error is sufficiently small %%%%%%%%%%%%%%%
 if(error(i)<1.5)
    cL_new=cL_old+rk4_RDM(cL_old,lD,kon,koff,dtc);
    delt(i)=dtc;
    time(i+1)=time(i)+dtc;
    i=i+1;
   
%%%%%% replace old with new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cL_old=cL_new;
    cL_arr(:,i)=cL_old;
%%%%%from profile of LGN calculate force on nucleus and centrosome
  
  if (xC_lin(i-1)>=xN_lin(i-1))
     fspr(i)=(xC_lin(i-1)-xN_lin(i-1)-eps); 
     %%%% spring force calculated from linear coordinate
     
%      [fact(i),Ip,In]=fcen_per(xC_per(i-1),Lm,x,dx,cL_new); 
%      [fact(i),Ip,In]=fcen_per_interp(xC_per(i-1),Lm,x,dx,dx_nuc,cL_new);
     [fact(i),Ip,In]=Fact_periodic(xC_per(i-1),Lm,x,dx,dx_nuc,cL_new);
     %%%% active force calculated from periodic coordinate
     
     xN_lin(i)=xN_lin(i-1)+dtc*fspr(i)/tauc; 
     xC_lin(i)=xC_lin(i-1)+dtc*((fact(i)/taua)-fspr(i)/tauc); 
   
     elseif (xC_lin(i-1)<xN_lin(i-1))
     fspr(i)=(xC_lin(i-1)-xN_lin(i-1)+eps);
%      [fact(i),Ip,In]=fcen_per(xC_per(i-1),Lm,x,dx,cL_new); 
%      [fact(i),Ip,In]=fcen_per_interp(xC_per(i-1),Lm,x,dx,dx_nuc,cL_new); 
     [fact(i),Ip,In]=Fact_periodic(xC_per(i-1),Lm,x,dx,dx_nuc,cL_new);
     
     xN_lin(i)=xN_lin(i-1)+dtc*fspr(i)/tauc;
     xC_lin(i)=xC_lin(i-1)+dtc*((fact(i)/taua)-fspr(i)/tauc); 
  end
  
  %%%%%%% update periodic coordinates from linear coordinates.
  fac_N=xN_lin(i)/L;   
  %%%% integer part of fac_N is the number of periodic image. 
 
  if fac_N>=0
   int_N=floor(fac_N); del_N=fac_N-int_N;
   xN_per(i)=del_N*L;
   
  elseif fac_N<0
   int_N=ceil(fac_N);del_N=fac_N-int_N;
   xN_per(i)=L*(1+del_N);
  
  end
  
  
  fac_C=xC_lin(i)/L; 
 
  if fac_C>=0
   int_C=floor(fac_C); del_C=fac_C-int_C;
   xC_per(i)=del_C*L;
   
  elseif fac_C<0
   int_C=ceil(fac_C); del_C=fac_C-int_C;
   xC_per(i)=L*(1+del_C);
      
  end
  
  
  %%%%%%%% plot and movie
  if (mod(i,5)==0)
  
    set(0,'DefaultAxesFontSize',14,'DefaultTextFontSize',14);
    
    figure(1)
    subplot(1,3,1)
    hold off
  
    plot(x,cL_new,'-','linewidth',1)
    hold on
    dsize=60;
    s1=scatter(xN_per(i),0); hold on; s1.Marker='o';     %%%% nucleus position
    set(s1,'SizeData',dsize,'MarkerFaceColor',[0 0.75 0],'MarkerEdgeColor',[0 0.75 0])

    s2=scatter(xC_per(i),0); hold on; s2.Marker='o';     %%%% centrosome position
    set(s2,'SizeData',dsize,'MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0.8 0 0])
  
  
  %%%%%% MT extent
  s3=scatter(x(Ip),zeros(1,size(Ip,2))); s3.Marker='.'; %% positive arm
  set(s3,'SizeData',dsize,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
    
  s4=scatter(x(In),zeros(1,size(In,2))); s4.Marker='.'; %% negative arm
  set(s4,'SizeData',dsize,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
    
    
  xlabel(['x in lab frame']); ylabel(['c_L(x)'])
  xlim([0 L]); 
  title([' time=',num2str(time(i))]); 

  axis square
  drawnow  
  
    
  %%%%%% Range of chromosomes
  subplot(1,3,2)
  hold off;
  plot(x,kon,'.r','Markersize',10);  hold on
  plot(x,koff,'.b','Markersize',10) 
  xlim([0 L]);  axis square; 
  xlabel(['x(\mum)']);  title([' time=',num2str(time(i))])
  drawnow  
    
  
  %%%% forces in time   
    subplot(1,3,3)
    hold off
    plot(time,fact/taua,'b','linewidth',1.5,'MarkerSize',15); hold on
    plot(time,fspr,'r','linewidth',1.5,'MarkerSize',15); 
    xlim([0 tmax]);  axis square; 
    xlabel(['t(min)']);  title([' time=',num2str(time(i))])
    drawnow  
    
   
    F = getframe(gcf);
    open(mov);
    writeVideo(mov,F);
  end     
 
end 
end
% close(fig)
 close(mov);


 
%%%%%%%% Solution at regular time points %%%
 t_reg=0:3:tmax;
 [Xo,To]=meshgrid(time(1:end),x);
 [Xn,Tn]=meshgrid(t_reg,x);
 
 %%%  c LGN
 cL_fin=interp2(Xo,To,cL_arr,Xn,Tn);
 
 %%% forces
 fa_fin=interp1(time,fact,t_reg);
 fs_fin=interp1(time,fspr,t_reg);
 
 %%% linear and periodic positions
 xNlin_fin=interp1(time,xN_lin,t_reg);
 xClin_fin=interp1(time,xC_lin,t_reg);
 xNper_fin=interp1(time,xN_per,t_reg);
 xCper_fin=interp1(time,xC_per,t_reg);
 
 Mat_N=zeros(size(xNlin_fin,2),2);
 Mat_C=zeros(size(xClin_fin,2),2);
 Mat_N(:,1)=xNlin_fin;  Mat_N(:,2)=xNper_fin;
 Mat_C(:,1)=xClin_fin;  Mat_C(:,2)=xCper_fin;
 
 dlmwrite(fullfile(dir,['cL.dat']),cL_fin,'delimiter','\t')
 dlmwrite(fullfile(dir,['f_act.dat']),fa_fin,'delimiter','\t')
 dlmwrite(fullfile(dir,['f_spr.dat']),fa_fin,'delimiter','\t')
 dlmwrite(fullfile(dir,['xN_pos.dat']),Mat_N,'delimiter','\t')
 dlmwrite(fullfile(dir,['xC_pos.dat']),Mat_C,'delimiter','\t')
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 






