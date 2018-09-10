%%%% load xgrid, cL_arr and phi_nuc_reg,phi_sp_reg for this. 
%%%% Also need cell_len and cell_wid to know the x position corresponding
%%%% to angles. 
plt_kym =1; 
opt_movie =1; 
load('phi_in=70_mu0=6eminus4.mat')
mov = VideoWriter(['Rounding_with_Ran_movie.avi']);

xint = xgrid(1,:)/xgrid(1,end); 
cnorm = zeros(size(cL_arr,1),size(cL_arr,2)); 
cL_arr = cL_arr*koffin;   %%% multiply by koffout to get density in units of kon/koffin.
for ifr =1:size(xgrid,1)
    %% define normalised x coordinate 
    xnorm = xgrid(ifr,:)/xgrid(ifr,end); 
    cnorm(:,ifr) = interp1(xnorm,cL_arr(:,ifr),xint); 
    
    ecc = sqrt(1-(cell_wid(ifr)^2)/cell_len(ifr)/cell_len(ifr));    
    shape = @(th) cell_wid(ifr)./sqrt(1-(ecc^2)*cos(th).^2);  %%%%  define ellipse
    th = round(linspace(0,2*pi,360*10),rn);
    dth = 2*pi/(size(th,2)-1);
    drdth = @(th) -cell_wid(ifr)*(ecc^2)*cos(th).*sin(th)./((1-(ecc^2)*(cos(th)).^2).^1.5);
    fun = @(th) sqrt(shape(th).^2 + drdth(th).^2); 
    Len_chk(ifr) = integral(fun,0,2*pi);
    xend(ifr) = xgrid(ifr,end);
    
    x_cen1(ifr) = integral(fun,0,phi_sp_reg(ifr))/Len_chk(ifr);
    x_nuc1(ifr) = integral(fun,0,phi_nuc_reg(ifr))/Len_chk(ifr);
    x_cen2(ifr) = integral(fun,0,phi_sp_reg(ifr)+pi)/Len_chk(ifr);
    x_nuc2(ifr) = integral(fun,0,phi_nuc_reg(ifr)+pi)/Len_chk(ifr);
    
    if opt_movie ==1
    tpnt = ifr;
    conc = cL_arr(:,ifr); 
    tim =  tarray(ifr); 
            
    cl_lo = 0; cl_hi = 5;
    if tim < Tround; Dynamic_plots_running_sim;  drawnow
    else; Dynamic_plots_running_sim_with_angle; drawnow
    end
            
    
        F = getframe(gcf);
        open(mov); 
        writeVideo(mov,F);
    end
    
    
end

if opt_movie ==1
    close(mov)
end

imagesc(tarray,xnorm,cnorm); colorbar
axis xy
hold on
plot(tarray,x_cen1,'.',tarray,x_cen2,'.'...
    ,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10)
plot(tarray,x_nuc1,'.',tarray,x_nuc2,'.'...
    ,'MarkerFaceColor','w','MarkerEdgeColor','w','MarkerSize',10)





