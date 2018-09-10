function [ft_form,ft_res,ft_gof]= AR_fits_DNA(quantity,norm_opt,fit_opt,col_mat,plt_indiv,plt_save)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% quantity = quantity to be plotted and fitted
%%%%%% norm_opt = option for normalisation, 1 for normalised, 0 for raw
%%%%%% fit_opt = option for fit : 1 for tanh function, 0 for exponential
%%%%%% col_mat = color matrix of size (Ncell, 3) specifyin RGB values for
%%%%%% each curve. plt_indiv=1 to make individual plots.
%%%%%% plt_save = 1 for saving plots, 0 otherwise.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ft_form =   functional form of the fit
%%%%%% ft_res = fit parameters
%%%%%% ft_gof = fit options


%%%%% data is extracted using the following code and the quantities are
%%%%% defined. 
Load_dynamic_data_All;

%%%% find the location of indices of new cell data
loc = [0 find(diff(Dynamic_shape_all.cell_id))' size(Dynamic_shape_all,1)]; 
% prompt = 'Enter the quantity to be plotted (fmt: Dynamic_shape_all.quantity): '; 
% %prompt = Dynamic_shape_all.dna_length./Dynamic_shape_all.dna_width; 
% 
% prompt2 = 'Enter 0 for raw plots, 1 for normalised plots: ';
quant_y = quantity; 
 
%%%% size of Lmat is set by the maximum duration data after NEB. The maximum
%%%% duration of the data in 58 cells is for a cell that lasts for 33
%%%% frames starting from NEB. 
Lmat = zeros(size(loc,2)-1,33);
%color_vec = 0.25*linspace(0,1,size(loc,2)-1);

R_vec = col_mat(:,1);
G_vec = col_mat(:,2);
B_vec = col_mat(:,3);


for icell = 1:size(loc,2)-1
    Lvec = zeros(1,size(Lmat,2));
    %%%% for check 
%   [Dynamic_shape.id(loc(icell)+1) Dynamic_shape.id(loc(icell+1)+1)]     
    time = Dynamic_shape_all.time(loc(icell)+1:loc(icell+1));    
    
    %%% find the 0 of time.
    t0 = find(time==0); 
    zero(icell) = t0; 
    quant_dyn = quant_y(loc(icell)+1:loc(icell+1));
    
    %%%% initial and final values of the quantity
    quant_I(icell) = quant_dyn(t0); 
    quant_F(icell) = quant_dyn(end); 
    tsize(icell) = size(time,1)-t0+1;
    
%     if plt_indiv ==1 
%         figure;
%     else 
%         figure(1)
%     end
    
    if norm_opt == 0 
        figure(1)
        
        plot(time(t0:end),quant_dyn(t0:end),'.-.','color',[R_vec(icell)  G_vec(icell) B_vec(icell)]...
            ,'linewidth',1.2,'Markersize',16)
        Lvec(1:tsize(icell)) = quant_dyn(t0:end); 
        Lmat(icell,:) = Lvec; 
        axis square
        set(gca,'fontsize',14)
        xlabel('Time','fontsize',16);
        hold on; drawnow

    elseif norm_opt ==1
        norm_y = (quant_dyn(t0:end)-quant_dyn(1))/(quant_F(icell)-quant_I(icell)); 
        plot(time(t0:end),norm_y,'.-.','color',[R_vec(icell)  G_vec(icell) B_vec(icell)],'linewidth',1.5)
        hold on; drawnow
        axis square
        set(gca,'fontsize',14)
        xlabel('Time','fontsize',16);

        Lvec(1:tsize(icell)) = norm_y(t0:end); 
        Lmat(icell,:) = Lvec; 
    
    end
    
    hold on
    if (plt_save==1 && plt_indiv ==1)      
        saveas(gca,['DNA_AR_Cell=',num2str(icell),'.png'])  
    elseif (plt_save==1 && plt_indiv ==1)
        saveas(gca,'DNA_AR_All_cells.png')
    end
    
    hold on;
    drawnow
%     close(gcf)
    
    hold on
    clear time 
    clear quant_dyn Lvec 
end

tvec = linspace(0,3*(max(tsize)-1),max(tsize));

for it = 1:size(tvec,2)
     I = find(Lmat(:,it)); 
     size_pnt(it) = size(I,1);
     mean_quant(it) = mean(Lmat(I,it));
     std_quant(it) = std(Lmat(I,it)); 
end

h1 = errorbar(tvec,mean_quant,std_quant,'ro','linewidth',2);
hold on 

in_frame =1; %%%%%% (NEB)
fin_frame = 14; %%%%% final frame upto which the fit is made. 

tnew = tvec(in_frame:fin_frame);
quant_new = mean_quant(in_frame:fin_frame); 
if fit_opt ==1   %%%%% hill function fit
     [xData, yData] = prepareCurveData( tnew, quant_new );
     % Set up fittype and options.
     ft_form = fittype( 'a*tanh((x-C)/tau)+b', 'independent', 'x', 'dependent', 'y' );
     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
     opts.Display = 'Off';
     opts.StartPoint = [0.45 13 10 1.95];

     % Fit model to data.
     [ft_res, ft_gof] = fit( xData, yData, ft_form, opts );
     
     ytanh = ft_res.a*tanh((tnew-ft_res.C)/ft_res.tau) + ft_res.b;
    
     hold on
     plot(tnew,ytanh,'k','linewidth',2)

else  %%%%% exp fit from t=6 mins
    [xData, yData] = prepareCurveData( tnew(3:end), quant_new(3:end) );

    %Set up fittype and options.
    ft_form = fittype( 'a*(1-exp(-b*x))+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.251083857976031 0.616044676146639 0.473288848902729];

    % Fit model to data.
    [fit_res, fit_gof] = fit( xData, yData, ft, opts );
    yexp = fit_res.a*(1-exp(-fit_res.b*tnew(3:end)))+ fit_res.c;

    hold on     
    plot(tnew(3:end),yexp,'linewidth',2)
    
end





%%%%%% equation for cell length 
%%%% cell_length(t) = 22.83 + q1.a*exp(q1.b*t);
%%%% [q1.a,q1.b] = [14.84, -0.24]

%%%%%% equation for cell width 
%%%% cell_width(t) = 22.35 + q1.a*exp(q1.b*t);
%%%% [q1.a,q1.b] = [1.62, -0.57]



























