%% load data 
%%%%% MAKE MOVIES AND PLOTS OF ALL THE QUANTITIES TO CLARIFY THE DOUBTS, 
%%% AND DEFINITIONS. 

Load_xls_file; 
id_all = analysis.cell_id; 
%% 
% %%  find locations of different cells
cell_loc = [0 find(diff(id_all))' size(id_all,1)];
plt_kym =1;
plt_save =1;
opt_write =1;
opt_movie =1;
opt_disp =1; 


% for icell = 1:size(cell_loc,2)-1
    icell  = 1
   
    if opt_movie ==1
    mov = VideoWriter(['Bipolar_Cell=',num2str(icell),'.avi']) ; 
    end
    %%%%%  extract data for each cell
    %%%%% Quantities constant at one time point
    %%% frame number
    frame_icell = analysis.frame(cell_loc(icell)+1: cell_loc(icell+1)); 
    %%% time points
    time_icell = analysis.time(cell_loc(icell)+1: cell_loc(icell+1)); 
    %%% Major axis
    major_icell = analysis.major(cell_loc(icell)+1:cell_loc(icell+1));
    %%% Cell length
    cell_length = analysis.length(cell_loc(icell)+1:cell_loc(icell+1));
    %%% minor axis
    minor_icell =  analysis.minor(cell_loc(icell)+1:cell_loc(icell+1));
    %%% mean lgn.
    mean_lgn_icell = analysis.lgn_mean(cell_loc(icell)+1:cell_loc(icell+1));
    %%%  z-height
    z_height = analysis.z(cell_loc(icell)+1:cell_loc(icell+1)); 
    %%% cell-angle
    cell_ang = analysis.initial_cell_angle(cell_loc(icell)+1:cell_loc(icell+1)); 
    %%% spindle angle
    sp_ang_icell = mod(analysis.theta(cell_loc(icell)+1:cell_loc(icell+1))+360,360); 
    
    %%%% Quantities that are s/angle dependent
    ang_coord = analysis.degrees(cell_loc(icell)+1:cell_loc(icell+1)); 
    s_coord = analysis.s(cell_loc(icell)+1:cell_loc(icell+1)); 
    %%%  Do not use s_norm_coord :  some issues with cut-off in the data.
    %%%  Use s_coord instead. 
    s_norm_coord = analysis.s_norm(cell_loc(icell)+1:cell_loc(icell+1));
    lgn_norm =  analysis.lgn_norm(cell_loc(icell)+1:cell_loc(icell+1));
    
    %%%% get the frame location :  determine starting points from this but
    %%%% end point for a fixed frame by the criteria : s < length
    frvec = find(diff(frame_icell));   
    fr_loc = [0 frvec' size(frame_icell,1)];
  
    %%%%% get all quantities at each frame 
    frame_vec = frame_icell(fr_loc(1:end-1)+1); 
    t_arr = time_icell(fr_loc(1:end-1)+1);
    major_ax_vec = major_icell(fr_loc(1:end-1)+1);
    minor_ax_vec = minor_icell(fr_loc(1:end-1)+1);
    cell_len_vec = cell_length(fr_loc(1:end-1)+1);
    mean_lgn_vec = mean_lgn_icell(fr_loc(1:end-1)+1); 
    z_vec = z_height(fr_loc(1:end-1)+1); 
       
    Time_table = table(frame_vec, t_arr, major_ax_vec,...
            minor_ax_vec, mean_lgn_vec, z_vec, cell_len_vec);
    
    %%%% To store lgn data 
    diff_vec = diff(fr_loc); 
    s_len = max(diff_vec); 
    
    %%%% initialise a matrix : lgn_mat = [s_coord ang_coord lgn_norm]
    %%%%  at each time point the length of these vectors is different,
    %%%%  replace the unavailable points with nan. Take care of this when 
    %%%%  reading this data file
    
    lgn_mat  = zeros(s_len,3*size(t_arr,1)); 
    lgn_mat(1:end,1:end) = nan; 
    sp_ang_vec = []; 
    for ifr =1:size(fr_loc,2)-1
    %%%% identify those points that are smaller than the cell length
%     Ivec = find(s_coord(fr_loc(ifr)+1:fr_loc(ifr+1)) <= cell_len_vec(ifr));  
    %%% This still did not work. Use the following criterion.
    
    %%%% index of the first point of the angular region thats counted twice
    %%%% since the angle coordinate changes from ~ 360 to something close
    %%%% 0, derivative at this point is negative. 
%     vec =  ang_coord(fr_loc(ifr)+1:fr_loc(ifr+1)); 
%     ind = find(diff(vec)<0); 
%     if size(ind,1) <1
%         ind = fr_loc(ifr+1) - fr_loc(ifr); 
%     end
    %%%% points to be deleted from the end since they are in twice measured
    %%%% region
    %miss_pnt = size(fr_loc(ifr)+1:fr_loc(ifr+1),2)-Ivec(end); 
    %%%% for check :  Ang_end should be ~= 360.
    Ang_int(ifr,1) = ang_coord(fr_loc(ifr+1));
    Ang_int(ifr,2) = ang_coord(fr_loc(ifr+1)); 
    s_int(ifr,1) = s_coord(fr_loc(ifr)+1); 
    s_int(ifr,2) = s_coord(fr_loc(ifr+1)); 
    
    lgn_mat(:,(ifr-1)*3+1) = s_coord(fr_loc(ifr)+1:fr_loc(ifr+1)); 
    lgn_mat(:,(ifr-1)*3+2) = ang_coord(fr_loc(ifr)+1:fr_loc(ifr+1)); 
    lgn_mat(:,(ifr-1)*3+3) = lgn_norm(fr_loc(ifr)+1:fr_loc(ifr+1)); 
   
    if t_arr(ifr) >= 0
        sp_ang_vec = [sp_ang_vec sp_ang_icell(fr_loc(ifr)+1)]; 
    end
    %%%%% plots and movies
        if opt_disp ==1
            dsize = 60;
            hold off
            plot(s_coord(fr_loc(ifr)+1:fr_loc(ifr)+ind), ...
                    lgn_norm(fr_loc(ifr)+1:fr_loc(ifr)+ind),...
                      '.-.','color',[1 0 0],'linewidth',2 ,'Markersize',15)
            hold on; 
            
            %%%% spindle angle
            if t_arr(ifr) >=0
            s2=scatter(sp_ang_vec(ifr-find(t_arr==0)+1),0); hold on; s1.Marker='o';  %%% wrt lab
            set(s2,'SizeData',dsize,'MarkerFaceColor',[0 0.75 0],'MarkerEdgeColor',[0 0.75 0])
            end
            
            xlabel('Angular coord','FontSize',16);  ylabel('Normalised LGN ','FontSize',16);
            set(gca,'fontsize',14)
            xlim([0 360]);  ylim([0 max(lgn_norm)+1])
            title(['Time =', num2str(t_arr(ifr))])
            if opt_movie ==1         
               F = getframe(gcf);
               open(mov);
               writeVideo(mov,F);  
            end
        end

    end
    if opt_movie ==1
        close(mov)
    end
    %%%%%%% Clarify data from Andrea
    
% %     if plt_kym ==1
% %        imagesc(t_arr,ang_arr,lgn_mat);  colorbar
% %        axis xy
% %        hold on 
% %        plot(t_arr,sp_ang_vec,'^','MarkerEdgeColor','k','MarkerFaceColor','k')
% %        xlabel('Time','FontSize',16);  ylabel('Linear coordinate','FontSize',16);
% %        set(gca,'FontSize',14)
% %        if plt_save ==1
% %           saveas(gcf,['Bipolar_rounding_Kym_Cell=',num2str(icell),'.png'])
% %        end
% %     end
    
    %%%%  write data for each cell
    if opt_write ==1
        dlmwrite(['Bipolar_LGN_kym_Cell=',num2str(icell),'.csv'],lgn_mat,'delimiter','\t')
        writetable(Time_table,['Bipolar_time_quants_Cell=',num2str(icell),'.dat'],'delimiter','\t')
    end
    
    
%  %%%%%  clear cell specific variables
%  clear lgn_mat Time_table ;
%  clear time_icell major_icell minor_icell sp_ang_icell mean_lgn_icell z_height cell_ang ;
%  clear ang_coord s_coord s_norm_coord lgn_norm frvec fr_loc s_size; 
%  clear sys_L Ang_end diff_vec s_max frame_vec t_arr minor_ax_vec major_ax_vec;
%  clear mean_lgn_vec z_vec sp_ang_vec
%  
%  
 
 
 
 
 
 
 
 
    
% end