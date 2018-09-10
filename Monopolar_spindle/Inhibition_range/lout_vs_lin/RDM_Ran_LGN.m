clear all; close all;
%%%%% This code checks the relation between Ran inhibition length scale and
%%%%% LGN rise length scale for different diffusion constants. The length
%%%%% of the system has been taken to correspond to Cell1. 
icell = 1; 
dir = (['Cell',num2str(icell)]); 
lgn_init = load(fullfile(dir,['Cell',num2str(icell),'_lgn_init.dat']));
xlab=lgn_init(:,1); 

%%%%%%% system and corrds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx=0.01; L=floor(xlab(end)/dx)*dx;  
x=(0:dx:L); N=size(x,2); 
tmax=200; dtc=0.01;                  %%% initial choice of dtc 

dtreg=1;  tol=0.001; 
dtm=dtc;                %%%%%% choice for the time step of  mechanical part

%%%%%%%%% Physical parameters of the system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
konout=1; konin=1;  koffout=0.02; koffin=0.1;         % in and out rates
dim=2;

% mov=VideoWriter(['D=',num2str(D),'va=',num2str(10*va),'.avi']);

%%%%%% Regular time array 
t_reg = 0:dtreg:tmax; 


Dvec = [0.001 0.01 0.1];
lvec = 5:15; 
jD = 3;
jl = 1;
% for jD = 1:size(Dvec,2)
% for jl = 1:size(lvec,2)
    l = lvec(jl); 
    D = Dvec(jD); 

i=1;  time(i)=0;   jreg=2;
delt(i) = dtc;
cL_old= rand(1,N);

%%%% Variables used for storing data at regular time points. Initialiase
%%%% for jreg=1.
cL_arr=zeros(N,size(t_reg,2));
cL_arr(:,1)=cL_old;

tic

if mod(N,2) ==0; xN = x(N/2) ;
else; xN = x((N+1)/2); end


%%%%%% update kon according to the nucleus position
kon=Rates(konin,konout,x,dx,xN,l,L); 
koff=Rates(koffin,koffout,x,dx,xN,l,L);
lD = sqrt(D./koff);
 
while (time(i)<tmax)   
    
    cL1=cL_old+rk4_RDM(cL_old,lD,kon,koff,dtc,dx);
    cLhalf=cL_old+rk4_RDM(cL_old,lD,kon,koff,0.5*dtc,dx);
    cL2=cLhalf+rk4_RDM(cLhalf,lD,kon,koff,0.5*dtc,dx);
   
    erru=norm(cL1-cL2);               
    error(i)=erru/tol; 
   
    if error(i)~=0
       fac=1/(error(i)^0.2);
       dtc=min([dtc*fac 10*dtc]);   %%%% time step update
    end
   
    %%%%% iterate solution if error is sufficiently small %%%%%%%%%%%%%%
    %%%%% and perform interpolation
    if(error(i)<1.5)
        cL_new=cL_old+rk4_RDM(cL_old,lD,kon,koff,dtc,dx);
        delt(i)=dtc;       
        time(i+1)=time(i)+dtc;
         
        %%%% check if t_reg(jreg) falls between time(i) and time(i-1)
        if (time(i)<t_reg(jreg) && t_reg(jreg)<=time(i+1))
        t_reg(jreg)
        Dif = abs(time(i+1)-time(i)); 
           
        if Dif >= 1e-3 
        t_irreg = [time(i) time(i+1)]; 

        ck_int(:,1) = cL_old; 
        ck_int(:,2) = cL_new; 

        [Xo,To] = meshgrid(t_irreg,x);
        [Xn,Tn] = meshgrid(t_reg(jreg),x); 
        
        cL_arr(:,jreg) = interp2(Xo,To,ck_int,Xn,Tn); 
        
        
        else
        cL_arr(:,jreg) = cL_new;      
         
        end
        
        jreg = jreg+1;
        end
           
%         plot(cL_new);  drawnow
        %%% advance i
        i=i+1;
        %%%%%% replace old with new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cL_old=cL_new;  
    end  
   
    
end

runtime = toc;
save(['D=',num2str(D),'l=',num2str(l),'.mat'])
% clear time  delt cL_new cL_old cL_arr

% 
% end 
% end 


 
 
 
 
 
 
 
 
 
 






