%%%%% This code plots Integral calculated from experimental
%%%%%  concentration vs simulated concentration for all cells. 

clear all;  close all ; 
sp_vec=[13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5];

option = 1;        %%%%% select distribution of MT ends : 1 for PDF, 2 for 1-CDF, 0 for step function.
opt_norm = 1;
plt_indiv = 0;
plt_save = 0;
dx = 0.02;        %%%%% space grid for LGN density
dx_nuc = dx;       %%%%% space grid for P(MT_ends). 
l = 11;            %%%%% inhibition range.

r = 2;
D = 0.01;
% param_data = load(['Optim_params_D=',num2str(D),'_tol=eminus3.dat']); 
% param_all = mean(param_data,1);

%%%%% Load probability density of MT ends
Mat = load('Prob_dist_MTends.dat'); 
MT_x = Mat(:,1);
MT_pdf = Mat(:,2);
param_all = load('history_D=0.01_niter=1_dx=02.dat');
Lm =20; 

fsim_all = [];
fexp_all = [];

for icell=1:size(sp_vec,2)
%     icell=1;
    sp_len = sp_vec(icell); 
    dir  = (['Cell',num2str(icell)]);     
    param = param_all(end,:);
    
    %%%% Call React_diff_solution to get tilde{c}. 
    [cL_arr,err] = React_diff_solution(dir,icell,sp_len,l,D,param,dx);    
    cL_arr = cL_arr*param(2);
    %%%%% experimental data  %%%%%%%%%%%%%%%%%%%%%
    %%%%% LGN concentration,  interpolate it on the grid x.
    lgn_init = load(fullfile(dir,['Cell',num2str(icell),'_lgn_init.dat']));  
    xlab = lgn_init(:,1); 
    %lgn_in = lgn_init(:,2); 
    %lgn_av = movmean(lgn_in,10);     %%%% smoothen data 

    %%%% this sets system size to be an integer multiple of numerical gridsize 
    L=floor(xlab(end)/dx)*dx; 
    x=(0:dx:L); N=size(x,2); 

    %%%%%%% LGN in space time, smoothening and interpolation on numerical grid
    lgn_data=load(fullfile(dir,['Cell',num2str(icell),'_conc.dat']));
    Nstep=size(lgn_data,1);                    %%%%% total number of time steps
    for ifr=1:Nstep
        lgn_aver(ifr,:)=movmean(lgn_data(ifr,:),10);
        lgn_exp(ifr,:)=interp1(xlab,lgn_aver(ifr,:),x);
        lgn_exp(ifr,:)=lgn_exp(ifr,:);
    end  
    
    %%%%% dna and centrosome positions
    dna_data = load(fullfile(dir,['Cell',num2str(icell),'_Mech.dat'])); 
    tvec = dna_data(:,1);  

    cen_pos = dna_data(:,2); 
    angle = dna_data(:,4);
    cen_pos_lin = dna_data(:,3);
    xNexp = cen_pos + sp_len*cos(pi/180*angle);
    xN_pos_lin = cen_pos_lin + sp_len*cos(pi/180*angle);
    
    %%%% to make DNA position periodic in lab frame.
    for ifr = 1:size(tvec,1)
        if (xNexp(ifr) < xlab(1))
            xNexp(ifr) = xlab(end)-abs(xNexp(ifr));
        elseif (xNexp(ifr) > xlab(end))
            xNexp(ifr) = xNexp(ifr)-xlab(end);
        end
    end
  
    %%%%%  mean velocity and relative velocity calculated from linear positions
    mean_pos = 0.5*(cen_pos_lin+xN_pos_lin); 
    rel_pos  = 0.5*(cen_pos_lin-xN_pos_lin); 

    meanvel_exp = diff(mean_pos)/3; 
    relvel_exp = diff(rel_pos)/3; 
    
    lgn_norm = tvec(end)*xlab(end)*lgn_exp/sum(sum(lgn_exp)); 
    cL_norm = tvec(end)*L*cL_arr/sum(sum(cL_arr)); 
   
    dt=3; 
    
    fact_sim=zeros(1,size(tvec,1));
    fact_exp=zeros(1,size(tvec,1)); 

    
    for ij=1:size(tvec,1) 
        %%%%% this block is to correct for the mistake in measurements. Sometimes 
        %%%%% position of centrosome is greater than the x-range on which LGN is
        %%%%% measured.  e.g.  cen_pos(46)=72.11 while xlab(end)=72.03 for cell=1
        if cen_pos(ij)>x(end)
            cen_pos(ij) = x(end)-dx; 
        elseif (cen_pos(ij)<0)
            cen_pos(ij) = 0; 
        end
        
        %%%%  calculate integral from normalised simulated and experimental
        %%%%  concentrations 
        [fact_sim(ij),Ipos,Ineg,xper,Pdfper] = Fact_periodic_Mtdist(...
                               cen_pos(ij),MT_x, MT_pdf,dx_nuc,x,dx,r,cL_arr(:,ij));
                           
        [fact_exp(ij),Ipos,Ineg,xper,Pdfper] = Fact_periodic_Mtdist(...
                              cen_pos(ij),MT_x, MT_pdf,dx_nuc,x,dx,r,lgn_norm(ij,:));       
    end
    
    if plt_indiv ==1
            figure;
        else
            figure(1)
    end  
        
    fsim_all = [fsim_all fact_sim];
    fexp_all = [fexp_all fact_exp];
    plot(fact_sim,fact_exp,'.','MarkerSize',20);  hold on; 
    set(gca,'Fontsize',14)
    ylabel('Integral of normalised experimental LGN','FontSize',16)
    xlabel('Integral of normalised simulated LGN','FontSize',16)
    grid on 
    box on 
    drawnow
    
    if (plt_save ==1 && plt_indiv ==1)
    saveas(gcf,['Exp_vs_Sim',num2str(icell),'.png'])
    end
    
    %%%%%  clear all cell specific variables
    clear lgn_aver lgn_exp lgn_data lgn_init lgn_in lgn_av L  x xlab lgn_norm
    clear L  x xlab cL_arr treg
    clear tvec dna_data cen_pos angle cen_pos_lin xNexp xN_pos_lin
    clear mean_vel rel_vel Mat
    clear fact_sim fact_exp
 
end
%legend('Cell 1','Cell 2','Cell 3','Cell 4','Cell 5','Cell 6','Cell 7','Cell 8','Cell 9','Cell 10', 'Cell 11','Cell 12' )

if (plt_save ==1 && plt_indiv ~=1)     
    saveas(gcf,'Exp_vs_Sim_integral_All.png')
end




% dlmwrite('Mean_conc.dat',total_mean)

%%%%%% plot of difference in concentration from exp and simulation as
%%%%%% function of velocity

% lgn_norm = tvec(end)*xlab(end)*lgn_exp/sum(sum(lgn_exp)); 
% cL_norm = tvec(end)*L*cL_arr/sum(sum(cL_arr));
% for ifr =2:size(tvec,1)    
%     DeltaC(ifr) =  sqrt(sum((cL_norm(:,ifr)-lgn_norm(ifr,:)').^2));
% end
% plot(meanvel_exp,DeltaC(2:end),'o')




