% %%%% load xgrid, cL_arr and phi_nuc_reg,phi_sp_reg for this. 
% %%%% Also need cell_len and cell_wid to know the x position corresponding
% %%%% to angles. 
rs = 10; rn =8;
plt_kym =0; 
opt_movie = 0 ; 
opt_x = 0; 
dir = 'Mobility_data'; 
mu0_vec = [1 0.1 0.01 0.002 0.001 0.0001 1e-5]; 
phi_vec = rad2deg(linspace(0,pi/2,19));

%%% Assign physical parameters.
ArmR = 9;    %%% half length of metaphase plate
RanR = 11;   %%% Range of RAN inhibiton
Lm = 15;     %%% Length of MTs
len_sp = 7;  %%% half length of spindle
D = 0.01;    %%% Diffusion constant
%mu0 = 0.00063; %%% mobility
Taud = 9;     %%% Delay between NEP and force application on spindle
konout = 1; konin = 1;          %%% On rates outside and inside inhibiton range
koffout = 0.039; koffin = 0.172;  %%% off rates outside and inside. 
dna_R = 13;        %%%%%% radius of DNA(9) + inhibition range(4) before t= Tround
Tround =9; 
T_total = 100; 

% %%%% load common parameters for all data
% com_data = load(fullfile(dir,['Time_quants_mu0=',num2str(10),'Phi_in=',num2str(rad2deg(0)),'.csv']));
% sp_angle = zeros(size(com_data,1),size(phi_vec,2)); 
% tarray = com_data(:,1);
% 
% for jmu =1:size(mu0_vec,2)
% %jmu = 4
%     for iphi =1:size(phi_vec,2)
%     
%     mu0 = mu0_vec(jmu);
%     phi_vec(iphi);
%     sp_angle_data = load(fullfile(dir,['Time_quants_mu0=',num2str(mu0),'Phi_in=',num2str(phi_vec(iphi)),'.csv']));
%     sp_angle(:,iphi) = sp_angle_data(:,4); 
%        
%     figure(jmu+1);
%     plot(tarray,rad2deg(sp_angle(:,iphi)),'LineWidth',1.5); hold on 
%     xlim([0,30]); set(gca,'fontsize',14); xlabel('Time (min)', 'Fontsize',16)
%     ylabel('Spindle angle in degrees', 'Fontsize',16)
%     title(['\mu_0=',num2str(mu0)])
%     drawnow
%     end
%     
%     yv1 = sp_angle(19,:); 
%     yv2 = sp_angle(61,:); 
%     
%     disp_54(jmu) = interp1(rad2deg(sp_angle(19,:)),rad2deg(sp_angle(61,:)),54.67); 
%     figure(1)
%     %%%% plot displacement at 30 mins
%     plot(rad2deg(sp_angle(1,:)),rad2deg(sp_angle(61,:)),'o-.','LineWidth'...
%          ,2,'MarkerSize',8,'MarkerFaceColor','w'); hold on; 
%     xlim([0,90]); set(gca,'fontsize',14); xlabel('Spindle angle at 9 mins', 'Fontsize',16)
%     ylabel('Spindle angle at 30 mins', 'Fontsize',16); ylim([0,90]);
%     drawnow
%   
% end

%% specify mobility and initial angle to plot kymograph 
mov = VideoWriter('TEst.avi'); 
% if plt_kym ==1 
     prompt1 = 'Enter mobility mu0 : ';
     prompt2 = 'Enter initial spindle angle in degrees : '; 
     mu0 = input(prompt1) ; 
     in_ang = input(prompt2);

    %%% load data 
    time_data = load(fullfile(dir,['Time_quants_mu0=',num2str(mu0),'Phi_in=',num2str(in_ang),'.csv']));
    xgrid = load(fullfile(dir,['Xgrid_mu0=',num2str(mu0),'Phi_in=',num2str(in_ang),'.csv'])); 
    cL_arr = load(fullfile(dir,['Conc_mu0=',num2str(mu0),'Phi_in=',num2str(in_ang),'.csv'])); 
    
    Nx = size(xgrid,1); 
    Nt = size(time_data,1);
    tarray = time_data(:,1);
    cell_len = time_data(:,2);
    cell_wid = time_data(:,3);
    phi_sp_reg = time_data(:,4);
    phi_nuc_reg = mod(phi_sp_reg + pi/2,2*pi);
    
    
    xint = xgrid(1,:)/xgrid(1,end); 
    cnorm = zeros(size(cL_arr,1),size(cL_arr,2)); 
    cL_arr = cL_arr*koffin;   %%% multiply by koffout to get density in units of kon/koffin.
    
    th_grid = linspace(0,2*pi,360);
    c_on_thgrid = zeros(size(th_grid,2),Nt);
       
    for ifr =1:size(xgrid,1)
        ifr 
        x = xgrid(ifr,:);   %%%  an assignment needed to make movie
        %%%%%% define normalised x coordinate 
        xnorm = xgrid(ifr,:)/xgrid(ifr,end); 
        cnorm(:,ifr) = interp1(xnorm,cL_arr(:,ifr),xint); 
        %plot(cnorm(:,ifr)); drawnow
        ecc = sqrt(1-(cell_wid(ifr)^2)/cell_len(ifr)/cell_len(ifr));    
        shape = @(th) cell_wid(ifr)./sqrt(1-(ecc^2)*cos(th).^2);  %%%%  define ellipse
        
        %dth = 2*pi/(size(th,2)-1);
        drdth = @(th) -cell_wid(ifr)*(ecc^2)*cos(th).*sin(th)./((1-(ecc^2)*(cos(th)).^2).^1.5);
        fun = @(th) sqrt(shape(th).^2 + drdth(th).^2);
        
        s_in(1) = 0;

        for ith =2:size(th_grid,2)
            s_in(ith) = round(integral(fun,0,th_grid(ith)),rs); 
        end
        
        
        Itest = find(s_in <= xgrid(ifr,end)); 
        if size(Itest,2)  == size(th_grid,2) 
            %%%% all points of s_in are  within the range
            sout = s_in;
            thout = th_grid;
        else
            sout = zeros(1,size(Itest,2)+1); 
            thout = zeros(1,size(Itest,2)+1); 
            sout(1:end-1) = s_in(Itest); 
            sout(end) = xgrid(ifr,end);
            thout(1:end-1) = th_grid(Itest); 
            thout(end) = th_grid(end); 
        end
              
        c_on_sgrid = interp1(xgrid(ifr,:),cL_arr(:,ifr),sout);
        %%% sout and thout have one-2-one mapping. 
        c_on_thgrid(:,ifr) = interp1(thout,c_on_sgrid,th_grid); 
        
        %%%% check nan in interpolated vector
        Iv = find(isnan(c_on_sgrid(:,ifr))); 
        if size(Iv,1) >=1 
            break 
        end
              
        L(ifr) = integral(fun,0,2*pi);
        xend(ifr) = xgrid(ifr,end);
       
        x_cen1(ifr) = integral(fun,0,phi_sp_reg(ifr))/L(ifr);
        x_nuc1(ifr) = integral(fun,0,phi_nuc_reg(ifr))/L(ifr);
        x_cen2(ifr) = integral(fun,0,phi_sp_reg(ifr)+pi)/L(ifr);
        x_nuc2(ifr) = integral(fun,0,phi_nuc_reg(ifr)+pi)/L(ifr);

        if opt_movie ==1
        tpnt = ifr;
        conc = cL_arr(:,ifr); 
        tim =  tarray(ifr); 

        cl_lo = 0; cl_hi = 5;
        if tim < Tround; Dynamic_plots_running_sim;  drawnow
        else; Dynamic_plots_running_sim_with_angle; drawnow
        end

        end 
        
        F = getframe(gcf);
        open(mov); 
        writeVideo(mov,F);

        clear sout s_in thout sout c_on_sgrid
    end

if opt_movie ==1
    close(mov)
end

if plt_kym ==1
if opt_x ==1   %%%% y-axis is normalised x
    imagesc(tarray,xnorm,cnorm); colorbar
    axis xy
    hold on
    plot(tarray,x_cen1,'.',tarray,x_cen2,'.'...
        ,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10)
    plot(tarray,x_nuc1,'.',tarray,x_nuc2,'.'...
        ,'MarkerFaceColor','w','MarkerEdgeColor','w','MarkerSize',10)

else          %%%% y-axis is theta
    
imagesc(tarray,rad2deg(th_grid),c_on_thgrid); colorbar
axis xy
hold on
plot(tarray,rad2deg(phi_sp_reg),'.',tarray,rad2deg(mod(phi_sp_reg+pi,2*pi)),'.'...
    ,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10)
plot(tarray,rad2deg(phi_nuc_reg),'.',tarray,rad2deg(mod(phi_nuc_reg+pi,2*pi)),'.'...
    ,'MarkerFaceColor','w','MarkerEdgeColor','w','MarkerSize',10)

    
end
end

%end


