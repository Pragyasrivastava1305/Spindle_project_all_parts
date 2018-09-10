%%%% code to produce the final figures as they appear in Fig.1. Run it 
%%%% sectionwise to create figures 1D

Load_dynamic_data_All; 
Ncell = 58; 

%% SECTION 1:  cell shape plots :  4 parameter fits
col_mat = zeros(Ncell,3); 
col_mat(:,1) = 0.5*ones(Ncell,1);
col_mat(:,2) = 1-0.25*linspace(0,1,Ncell);
figure(1)
[len_fit,len_opt,fig_hand] = Shape_plots_fits(Dynamic_shape_all.cell_length...
                        ,0,col_mat,0);
                    
hold on;
col_mat(:,1) = linspace(0.4,0.6,Ncell);
col_mat(:,2) = linspace(0.6,0.85,Ncell);
col_mat(:,3) = 0.9*ones(Ncell,1);
[wid_fit,wid_opt,fig_hand2] = Shape_plots_fits(Dynamic_shape_all.cell_width...
                        ,0,col_mat,0);
                    
                    
                    
                    
%% SECTION 2:  Cell length and width : 3 parameter tanh fits :  Fig. 1 D
%%%%% the range of variation of each other. The constant parameter in the
%%%%% fit function, in Param2_tanh_fit needs to be set to 22.6250 before
%%%%% calling it. 

col_mat(:,1) = 0.37*ones(Ncell,1); 
col_mat(:,2) = 0.79*ones(Ncell,1);
col_mat(:,3) = 0.75*ones(Ncell,1);

[len3_fit,len3_opt,fig_hand] = Param3_tanh_fit(Dynamic_shape_all.cell_length...
                        ,22.6250,0,col_mat,0);

col_mat(:,1) = ones(Ncell,1);
col_mat(:,2) = 0.8*ones(Ncell,1);
col_mat(:,3) = 0.64*ones(Ncell,1);

[wid3_fit,wid3_opt,fig_hand] = Param3_tanh_fit(Dynamic_shape_all.cell_width...
                        ,22.6250,0,col_mat,0);

                    
                    
%% SECTION 3: Cell AR plot: Fit to AR is forced to be 1 at large times. In reality,
%%%%% it is 1.05. The measured width and lengths are within the range of
%%%%% variation hence this is correct.  Produces Fig 1E. The constant parameter in the
%%%%% fit function, in Param2_tanh_fit needs to be set to 1 before
%%%%% calling it. 
col_mat = zeros(Ncell,3); 
col_mat(:,1) = 0.5*ones(Ncell,1);
col_mat(:,2) = 1-0.25*linspace(0,1,Ncell);
[AR_fit,AR_opt,fig_hand3]= Param3_tanh_fit(Dynamic_shape_all.cell_length./Dynamic_shape_all.cell_width...
                            ,1,0,col_mat,0); 

%% SECTION 4:  dna length and width plots
figure(2)
col_mat(:,1) = 0.37*ones(Ncell,1); 
col_mat(:,2) = 0.79*ones(Ncell,1);
col_mat(:,3) = 0.75*ones(Ncell,1);
[len_fit,len_opt,fig_hand] = Shape_plots_fits(Dynamic_shape_all.dna_length...
                        ,0,col_mat,0);
                    
hold on;
col_mat(:,1) = ones(Ncell,1);
col_mat(:,2) = 0.8*ones(Ncell,1);
col_mat(:,3) = 0.64*ones(Ncell,1);
[wid_fit,wid_opt,fig_hand2] = Shape_plots_fits(Dynamic_shape_all.dna_width...
                        ,0,col_mat,0);
 
                    
                    
%% DNA AR plots   :  4 parameter tanh fit. :  fit details in AR_fits_DNA, outputted in 
%%%%% ft_form
figure(3)                    
col_mat = zeros(Ncell,3); 
col_mat(:,1) = 0.5*ones(Ncell,1);
col_mat(:,2) = 1-0.25*linspace(0,1,Ncell);
DNA_AR = Dynamic_shape_all.dna_length./Dynamic_shape_all.dna_width; 

[ft_form,ft_res,ft_gof] = AR_fits_DNA(DNA_AR,0,1,col_mat,0,0); 
ylabel(['DNA Aspect Ratio'],'Fontsize',16)
box on; grid on; 























