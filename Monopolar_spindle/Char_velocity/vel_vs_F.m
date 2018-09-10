%%%% This code needs the data from Mono_F_vs_V_master_code.m, which writes
%%%% the data files being used in this code. Before running this code
%%%% always check if the data file being used is with the correct
%%%% normalisation of concentration. As taken into account in
%%%% Mono_F_vs_V_master_code.m, there are 4 options that can be used : 
%%%% (raw, exp), (raw, sim), (normalised, exp), (normalised,sim).
close all

plt_save =1;
plt_indiv = 0;
opt_exp = 0;    %%%%  1 for experimental data
opt_norm = 0;   %%%% 1 if normalised data

vel_vec_all = []; 
int_vec_all = [];

for icell =1:12
%     icell = 9
    dir  = (['Cell',num2str(icell)]);    
    %%%%% load appropriate file as specified by options in the beginning of
    %%%%% this code
    if (opt_exp ==1 && opt_norm ~=1)
        cdata = load(fullfile(dir,['Exp_raw_Force_velocity_PDF_cell=',num2str(icell),'.dat']));           
    elseif (opt_exp ==1 && opt_norm ==1)
        cdata = load(fullfile(dir,['Exp_norm_Force_velocity_PDF_cell=',num2str(icell),'.dat'])) ;
    elseif (opt_exp ~=1 && opt_norm ~=1)
        cdata = load(fullfile(dir,['Sim_raw_Force_velocity_PDF_cell=',num2str(icell),'.dat']));
    elseif (opt_exp ~=1 && opt_norm ==1)
        cdata = load(fullfile(dir,['Sim_norm_Force_velocity_PDF_cell=',num2str(icell),'.dat'])); 
    end
    
    %data1 = load(fullfile(dir,['Exp_Force_velocity_PDF_cell=',num2str(icell),'.dat']));
    vel_exp=cdata(:,1);     
    fact = cdata(:,2); 
    
    %%%%% make a single vector of xvalues and yvales to be used for fits
    vel_vec_all = [vel_vec_all; vel_exp]; 
    int_vec_all = [int_vec_all; fact];
   
    %%%% plot the data separately for each cell if plt_indiv =1
    if plt_indiv ==1
        figure;
    else
        figure(1)
    end
     
    plot(fact,vel_exp,'.','markersize',20,'linewidth',1);  hold on ;
    set(gca,'Fontsize',14)  
    ylabel('Observed velocity','FontSize',16)
    
    %%%%  title of plots
    if (opt_exp ==1 && opt_norm ~=1)       
        title('Integral of raw experimental LGN','FontSize',16)
%     axis square
    elseif (opt_exp ==1 && opt_norm ==1)
        title('Integral of normalised experimental LGN','FontSize',16)
    elseif (opt_exp ~=1 && opt_norm ~=1)
        title('Integral of raw simulated LGN','FontSize',16)
    elseif (opt_exp ~=1 && opt_norm ==1)
        title('Integral of normalised simulated LGN','FontSize',16)
    end
    
    grid on 
    box on
    drawnow
     
    %%%% save plots
    if (plt_save ==1 && plt_indiv ==1)    
    if (opt_exp ==1 && opt_norm ~=1)       
        saveas(gcf,['Exp_raw_Cell=',num2str(icell),'.png'])
    elseif (opt_exp ==1 && opt_norm ==1)
        saveas(gcf,['Exp_norm_Cell=',num2str(icell),'.png'])
    elseif (opt_exp ~=1 && opt_norm ~=1)
        saveas(gcf,['Sim_raw_Cell=',num2str(icell),'.png'])
    elseif (opt_exp ~=1 && opt_norm ==1)
        saveas(gcf,['Sim_norm_Cell=',num2str(icell),'.png'])
    end    
    end
   clear meanvel_exp fact
    
end 

hold on; 
%%%%%% plots of all cells
if (plt_save ==1 && plt_indiv ~=1) 
    
%%%%%  fits % Fit: 'linear fit' with no intercept
[xData, yData] = prepareCurveData( int_vec_all, vel_vec_all );
% Set up fittype and options.
ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.144008605579303;

% Fit model to data.
[fit_data, gof] = fit( xData, yData, ft, opts );
slope =fit_data.a; 
rsq = gof.rsquare; 
plot(int_vec_all,slope*int_vec_all,'k', 'linewidth',2) 
 
legend('Cell 1','Cell 2','Cell 3','Cell 4',...
        'Cell 5','Cell 6','Cell 7','Cell 8','Cell 9'...
        ,'Cell 10', 'Cell 11','Cell 12')
    if (opt_exp ==1 && opt_norm ~=1)       
        saveas(gcf,['Exp_raw_All.fig'])
    elseif (opt_exp ==1 && opt_norm ==1)
        saveas(gcf,['Exp_norm_All.fig'])
    elseif (opt_exp ~=1 && opt_norm ~=1)
        saveas(gcf,['Sim_raw_All.fig'])
    elseif (opt_exp ~=1 && opt_norm ==1)
        saveas(gcf,['Sim_norm_All.fig'])
    end
end
 







 
 

 
 
 
 
 
 
 
 
 
 
 
