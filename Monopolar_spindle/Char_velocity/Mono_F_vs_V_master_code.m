
%%%%% This program selects initial condition from each experiment to solve
%%%%% numerically for c_LGN, which is then used to calculate
%%%%% the integral of G(s-s')c_LGN(s') over s'. This integral is plotted
%%%%% against the obesrved velocity in experiments to extract
%%%%% characteristic velocity. The parameters of reaction diffusion
%%%%% equation are fixed at the optimised values. To make movies for
%%%%% individual cells, comment the for loop and specify icell.
%%%%% This code plots and writes several quantities. So it is possible to
%%%%% just manipulate this code to plot various quantities. But to make
%%%%% plots of (i) simulation vs experimental integral, (ii) integral vs
%%%%% obs velocity and (iii) total LGN in time from simulation and
%%%%% experiments, there are other small codes, which can be used. 

clear all;  close all ; 
sp_vec=[13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5];
option = 1;        %%%%% select distribution of MT ends : 1 for PDF, 2 for 1-CDF, 0 for step function.
opt_mov = 1;       %%%% set to 1 for movies
opt_exp = 0;       %%%% set to 1 if experimental LGN is to be used for movies and file writing
opt_write = 0;     %%%% set to 1 to write in files.
opt_norm = 1;      %%%% set to 1 if lgn concentration is to be normalised by its spatiotemporal mean
plt_ctime = 0;     %%%% set to 1 to generate plot of total lgn in time. But there is a separate code for this as well.
dx = 0.02;         %%%%% space grid for LGN density
dx_nuc = dx;       %%%%% space grid for P(MT_ends). 
l = 11;            %%%%% inhibition range.

r = 2;
D = 0.01;
param_data = load(['history_D=',num2str(D),'_niter=1_dx=02.dat']); 
minparam = param_data(end,:);

%%%%% Load probability density of MT ends
Mat = load('Prob_dist_MTends.dat'); 
MT_x = Mat(:,1);
MT_pdf = Mat(:,2);
Lm =20;   %%%%% length of MTs to be used if step distribution is chosen

%%%%% Load centrosome velocity from the new file by Andrea
load_centrosome_vel;
%%%% above code load cell_id and centrosome velocity in the two vectors:
%%%% cell_id, v_centrosome
id_loc = find(diff(cell_id));   %%%% pick up the indices where new cell_id starts.
id_loc = [0; id_loc; size(cell_id,1)];


for icell=1:size(sp_vec,2)
    icell 
    
    sp_len = sp_vec(icell); 
    dir  = (['Cell',num2str(icell)]);     
    param = [minparam(1:3),minparam(3+icell)];
    
    %%%% Call React_diff_solution to get \bar{c}, which is in unit of kon/koffnear.  
    [cL_arr,err] = React_diff_solution(dir,icell,sp_len,l,D,param,dx); 
    cL_arr = minparam(2)*cL_arr;        
    %%%%% experimental data  %%%%%%%%%%%%%%%%%%%%%
    %%%%% LGN concentration,  interpolate it on the grid x. 
    %%%%% lgn_init = [x_exp, lgn_initial(x_exp)], raw data from Andrea.
    
    lgn_init = load(fullfile(dir,['Cell',num2str(icell),'_lgn_init.dat']));  
    xlab = lgn_init(:,1);  
    vec = diff(xlab);
    dxlab =  vec(1);
    
    %%%% To smoothen data :  commented in most recent results, so no
    %%%% smoothening
%     lgn_in = lgn_init(:,2); 
%     lgn_av = movmean(lgn_in,10);     
    
    %%%% this sets system size to be an integer multiple of numerical gridsize 
    sys_L(icell) = xlab(end);   %%%% record experimental system size to later correlate with the max velocity
    L=floor(xlab(end)/dx)*dx; 
    x=(0:dx:L); N=size(x,2); 
    
    %%%%%%% LGN in space time, smoothening and interpolation on numerical grid
    %%%%% lgn_data is a space time matrix, with data on xlab space-grids,
    %%%%% and tvec (to be loaded from another file 'Cell_Mech.dat') time points.
    
    lgn_data=load(fullfile(dir,['Cell',num2str(icell),'_conc.dat']));
    Nstep=size(lgn_data,1);                    %%%%% total number of time steps
    
    for ifr=1:Nstep
%         lgn_aver(ifr,:)=movmean(lgn_data(ifr,:),10);
        lgn_exp(ifr,:)=interp1(xlab,lgn_data(ifr,:),x);
    end  
  
    %%%%% dna and centrosome positions
    dna_data = load(fullfile(dir,['Cell',num2str(icell),'_Mech.dat'])); 
    tvec = dna_data(:,1);     %%%%%  experimental time vector

    cen_pos = dna_data(:,2); 
    angle = dna_data(:,4);
    cen_pos_lin = dna_data(:,3);
    
    %%%% Calculate position of DNA (from centrosome position and angle). 
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
  
%   %%%%%  mean velocity and relative velocity calculated from PROJECTED DNA positions
    mean_pos = 0.5*(cen_pos_lin+xN_pos_lin); 
    rel_pos  = 0.5*(cen_pos_lin-xN_pos_lin); 
    meanvel_exp = diff(mean_pos)/3; 
    relvel_exp = diff(rel_pos)/3; 

    %%%%% assign correct velocities of centrosome, AS MEASURED BY ANDREA
    strt_ind = id_loc(icell)+1; 
    fin_ind = id_loc(icell+1);
    vel_cen =  v_centrosome(strt_ind+1:fin_ind).*sign(meanvel_exp);
    vmax(icell) = max(abs(vel_cen));
    vmin(icell) = min(abs(vel_cen));
    
    %%%% scatter plot of instantaneous velocity ranges in each cell.
    %%% scatter(icell*ones(size(vel_cen,1),1),abs(vel_cen)); hold on; drawnow
    
    %%%%% Total and mean spatiotemporal concentration. 
    exp_tot(icell) = sum(sum(lgn_exp))*dx*3; 
    exp_mean(icell) = exp_tot(icell)/tvec(end)/xlab(end);  %%%% grid spacing for lgn_exp is dx
    sim_tot(icell) = sum(sum(cL_arr))*dx*3; 
    sim_mean(icell) =  sim_tot(icell)/L/tvec(end);
    
    %%%% Total concentration integrated in space as a function of time 
    expC_as_t = sum(lgn_exp,2)*dx;
    simC_as_t = sum(cL_arr,1)*dx; 
    
    if plt_ctime ==1
        figure(1)
        plot(tvec,expC_as_t,'linewidth',1.5);  hold on
        figure(2);
        plot(tvec, simC_as_t,'linewidth',1.5); hold on 
    end
    
    lgn_norm = lgn_exp/exp_mean(icell); 
    cL_norm = cL_arr/sim_mean(icell); 
   
    %%%%%%% PART 2 :  calulation of forces from given lgn_norm and MT PDF for
    %%%%%%% each cell. 
    %%%%% initialise movie
    if option == 0
        mov = VideoWriter(['Force_Step_cell=',num2str(icell),'.avi']); 
    elseif option == 1
%         mov = VideoWriter(['Force_PDF_cell=', num2str(icell),'.avi']);
          if (opt_exp==1 && opt_norm ==1)
          mov = VideoWriter(['Cell=', num2str(icell),'PDF_Exp_Norm.avi']);
          elseif (opt_exp==1 && opt_norm ~=1) 
          mov = VideoWriter(['Cell=', num2str(icell),'PDF_Exp_Raw.avi']);  
          elseif (opt_exp ~=1 && opt_norm ==1)
          mov = VideoWriter(['Cell=', num2str(icell),'PDF_Sim_Norm.avi']); 
          elseif (opt_exp ~=1 && opt_norm ~=1)
          mov = VideoWriter(['Cell=', num2str(icell),'PDF_Sim_Raw.avi']);
          end
    elseif option ==2 
        MT_pdf = 1- cumsum(MT_pdf);
        if (opt_exp==1 && opt_norm ==1)
        mov = VideoWriter(['Cell=', num2str(icell),'CDF_Exp_Norm.avi']);
        elseif (opt_exp==1 && opt_norm ~=1) 
        mov = VideoWriter(['Cell=', num2str(icell),'CDF_Exp_Raw.avi']);  
        elseif (opt_exp ~=1 && opt_norm ==1)
        mov = VideoWriter(['Cell=', num2str(icell),'CDF_Sim_Norm.avi']); 
        elseif (opt_exp ~=1 && opt_norm ~=1)
        mov = VideoWriter(['Cell=', num2str(icell),'CDF_Sim_Raw.avi']);
        end      
    end

    fig=figure;
    epsilon=3.25;            %%%% separation between centrosome and nucleus
                             %%% calculated as mean distance between
                             %%% centrosome and DNA from cell 1. 
    dt=3; 
    fspr=zeros(1,size(tvec,1));
    fact=zeros(1,size(tvec,1));

    for ij=1:size(tvec,1) 
        
        %%%% assign concentration to be used to calculate integral. options
        %%%% are experimental (opt_exp==1), normalised by spatiotemporal mean 
        %%%% (opt_norm==1), simulation (opt_exp~=1) and raw (opt_norm~=1) 
        if (opt_exp ==1 && opt_norm ==1)
           cL_new = lgn_norm(ij,:);            
        elseif (opt_exp ==1 && opt_norm ==0)
           cL_new=lgn_exp(ij,:);     
        elseif (opt_exp ~=1 && opt_norm ==1)
           cL_new=cL_norm(:,ij);           
        elseif (opt_exp ~=1 && opt_norm ~=1)
           cL_new=cL_arr(:,ij);     
        end
                   
        %%%%% Use linear positions to calculate spring force 
        if (cen_pos_lin(ij) >= xN_pos_lin(ij))
           fspr(ij)=(cen_pos_lin(ij)-xN_pos_lin(ij)-eps); 
        elseif (cen_pos_lin(ij) < xN_pos_lin(ij))
           fspr(ij)=(cen_pos_lin(ij)-xN_pos_lin(ij)+eps);
        end  
     
 %%%%% this block is to correct for the mistake in measurements. Sometimes 
 %%%%% position of centrosome is greater than the x-range on which LGN is
 %%%%% measured.  e.g.  cen_pos(46)=72.11 while xlab(end)=72.03 for cell=1
         if cen_pos(ij)>x(end)
            cen_pos(ij) = x(end)-dx; 
        elseif (cen_pos(ij)<0)
            cen_pos(ij) = 0; 
        end
  
 %%%% active force calculated from periodic coordinate
        if option ~=0
        [fact(ij),Ipos,Ineg,xper,Pdfper] = Fact_periodic_Mtdist(...
                                cen_pos(ij),MT_x, MT_pdf,dx_nuc,x,dx,r,cL_new);
        else
        [fact(ij),Ipos,Ineg] = Fact_periodic_MTstep(cen_pos(ij),Lm,x,dx,dx_nuc,cL_new); 
        end
   
        %[fact_sim(icell,ij),Ipos,Ineg,xper,Pdfper] = Fact_periodic_Mtdist(...
        %                       cen_pos(ij),MT_x, MT_pdf,dx_nuc,x,dx,r,cL_norm(:,ij));
        dsize=60;
        
        if (option ~=0 && opt_mov ==1)    %%%% make movie when MT distribution is given from experiments.
        %subplot(1,2,1)
  
        hold off
      
        %%%% MT distribution in periodic coordinates
        s0=scatter(xper,Pdfper); hold on;  s0.Marker ='o';
        set(s0,'SizeData',0.2*dsize,'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[0 0 0])
  
        s1=scatter(cen_pos(ij),0) ;   s1.Marker = 'o'; 
        set(s1,'SizeData',dsize,'MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0 0 0])
  
        s2=scatter(xNexp(ij),0) ;   s2.Marker = 'o'; 
        set(s2,'SizeData',dsize,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0 0])
  
        hold on 

        plot(x,cL_new, '-g', 'linewidth',1.5);  hold on
        %%%%% positive branch
        s2=scatter(x(Ipos),cL_new(Ipos),'linewidth',2);   s2.Marker = '.'; 
        set(s2,'SizeData',dsize,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
        hold on 
         %    drawnow
        %%%%% negative branch
        s3=scatter(x(Ineg),cL_new(Ineg),'linewidth',2);   s3.Marker = '.'; 
        set(s3,'SizeData',dsize,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0])
  
        box on;  grid on
        hold on
        
        legend('MT pdf','Centrosome','DNA', 'LGN profile','LGN in pos region', 'LGN in neg region')
        xlim([0 xlab(end)]);    ylim([0 0.02+max(max(lgn_norm))]); 
        title(['Time=',num2str(tvec(ij)),'    Ratio of max to mean='...
                ,num2str(max(max(lgn_exp))/exp_mean(icell))])
        xlabel('Distance along cell periphery (\mum)')
        ylabel('LGN density')
        
        drawnow
  
        Fprime = getframe(fig);
        open(mov);
        writeVideo(mov,Fprime);  
   
        end
%         
%         if ij ==2 
%             stop
%         end 
    end
    
    if option ~= 0
    close(mov) 
    close(fig)
    end

   % Mat = [meanvel_exp relvel_exp fact(1:end-1)' fspr(1:end-1)'];
    Mat = [vel_cen fact(1:end-1)' fspr(1:end-1)'];

    %%%%%% write data files 
    if opt_write ==1
    if option == 0
        if (opt_exp ==1 && opt_norm ~=1)
            dlmwrite(fullfile(dir,['Exp_raw_Force_velocity_step_cell=',num2str(icell),'.dat']),Mat)           
        elseif (opt_exp ==1 && opt_norm ==1)
            dlmwrite(fullfile(dir,['Exp_norm_Force_velocity_step_cell=',num2str(icell),'.dat']),Mat) 
        elseif (opt_exp ~=1 && opt_norm ~=1)
            dlmwrite(fullfile(dir,['Sim_raw_Force_velocity_step_cell=',num2str(icell),'.dat']),Mat)
        elseif (opt_exp ~=1 && opt_norm ==1)
            dlmwrite(fullfile(dir,['Sim_norm_Force_velocity_step_cell=',num2str(icell),'.dat']),Mat) 
        end
        
    elseif option == 1 
        if (opt_exp ==1 && opt_norm ~=1)
            dlmwrite(fullfile(dir,['Exp_raw_Force_velocity_PDF_cell=',num2str(icell),'.dat']),Mat)           
        elseif (opt_exp ==1 && opt_norm ==1)
            dlmwrite(fullfile(dir,['Exp_norm_Force_velocity_PDF_cell=',num2str(icell),'.dat']),Mat) 
        elseif (opt_exp ~=1 && opt_norm ~=1)
            dlmwrite(fullfile(dir,['Sim_raw_Force_velocity_PDF_cell=',num2str(icell),'.dat']),Mat)
        elseif (opt_exp ~=1 && opt_norm ==1)
            dlmwrite(fullfile(dir,['Sim_norm_Force_velocity_PDF_cell=',num2str(icell),'.dat']),Mat) 
        end
       
    elseif option ==2    
        if (opt_exp ==1 && opt_norm ~=1)
            dlmwrite(fullfile(dir,['Exp_raw_Force_velocity_CDF_cell=',num2str(icell),'.dat']),Mat)           
        elseif (opt_exp ==1 && opt_norm ==1)
            dlmwrite(fullfile(dir,['Exp_norm_Force_velocity_CDF_cell=',num2str(icell),'.dat']),Mat) 
        elseif (opt_exp ~=1 && opt_norm ~=1)
            dlmwrite(fullfile(dir,['Sim_raw_Force_velocity_cDF_cell=',num2str(icell),'.dat']),Mat)
        elseif (opt_exp ~=1 && opt_norm ==1)
            dlmwrite(fullfile(dir,['Sim_norm_Force_velocity_CDF_cell=',num2str(icell),'.dat']),Mat) 
        end
               
    end
    end
             
    %%%%%  clear all cell specific variables
    clear lgn_aver lgn_exp lgn_data lgn_init lgn_in lgn_av L  x xlab lgn_norm
    clear L  x xlab cL_arr treg
    clear tvec dna_data cen_pos angle cen_pos_lin xNexp xN_pos_lin
    clear mean_vel rel_vel Mat
 
end

Mat_conc = [ exp_tot' exp_mean' sim_tot' sim_mean']; 
dlmwrite('Conc_mean_n_total.dat',Mat_conc)

%%%%%% plot of difference in concentration from exp and simulation as
%%%%%% function of velocity

% lgn_norm = tvec(end)*xlab(end)*lgn_exp/sum(sum(lgn_exp)); 
% cL_norm = tvec(end)*L*cL_arr/sum(sum(cL_arr));
% for ifr =2:size(tvec,1)    
%     DeltaC(ifr) =  sqrt(sum((cL_norm(:,ifr)-lgn_norm(ifr,:)').^2));
% end
% plot(meanvel_exp,DeltaC(2:end),'o')




