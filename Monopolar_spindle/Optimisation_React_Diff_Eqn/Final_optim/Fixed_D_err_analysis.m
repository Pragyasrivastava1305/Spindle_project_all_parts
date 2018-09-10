%%%%%%%  This code has three sections 
%%% First section loads the data on optimised parameter values for all
%%% iterations (20) for D = [1]. It then calls
%%% Error_all_cells for each row of parameters and creates a matrix err
%%% with 4 rows and 20 columns. If plt_err =1 then the section creates the
%%% figure of Error as function of different iterations and for all 4 Ds. 

%%% Second section plots parameters from all iterations for D =0.01 (since
%%% error minimum appears for D =0.01) and saves the figure. 

%%% Last section plots kymographs for parameters obtained from the mean of
%%% all equivalent optimisation results.

ncount =0; 
Dvec = [0 0.01 0.1 1];
param_small = []; 
Tot_iter = 20; 
plt_err =0;  %%%% 1 for error as function of iteration plot
plt_min =0;  %%%% 1 for parameter scatter plot for min D. 
plt_rates = 0; %%%% 1 if 3 parameter scatter plot
plt_all = 0;   %%%% 1 for all parameters scatter plot (rates + alpha)
plt_kym = 1;  %%%% 1 for all cell kymographs if 
plt_save =1;  %%%% 1 for saving the plot
dxvec = [0.05 0.05 0.1 0.5]; 
l = 11;

%for iD =1:size(Dvec,2)
iD = 4; 
D = Dvec(iD);
for iter=1:Tot_iter
        [iD iter]
        data = load(['history_D=1_niter=',num2str(iter),'.dat']); 
        param= data(end,:);
        err(iter) = Error_all_cell_Jun18(l,D,param,dxvec(iD));        
end
    
if plt_err ==1
plot(1:Tot_iter,err,'o','linewidth',1.5); hold on ; drawnow
end
%end

axis square; 
grid on; box on;
set(gca,'linewidth',1,'Fontsize',14)
set(gca,'yscale','log')
xlabel('Iteration number','Fontsize',16)
ylabel('Error','Fontsize',16) 
%legend('D=0','D=0.01','D=0.1','D=1')

if (plt_save ==1 && plt_err ==1)
    saveas(gcf,'Error_all_D.fig')
    saveas(gcf,'Error_all_D.pdf')    
end


%%
rc = linspace(0,1,20); 
figure; clc

if plt_min ==1
%%%%%  minimum appears for D=1, iter=1;
%%%%% generate final kymographs 
D = 1;
data = load(['Optim_params_D=',num2str(D),'_tol=eminus3.dat']); 
for iter =1:20

if plt_rates ==1
    plot(1:3,data(iter,1:3),'o','color',[rc(iter) 0 0],'Markersize',8,'linewidth',2); drawnow;  hold on   
    elseif plt_all ==1
    plot(1:15,data(iter,1:end),'o','color',[rc(iter) 0 0],'Markersize',8,'linewidth',2); drawnow;  hold on
end

end
axis square
if plt_rates ==1
    xticks(1:3)
    xticklabels({ 'k^{off}_{out}', 'k^{off}_{in}','k^{on}_{out}'})
    elseif plt_all ==1
    xticks(1:15)
    xticklabels({'k^{off}_{out}', 'k^{off}_{in}','k^{on}_{out}',...
        '\alpha_1' ,  '\alpha_2', '\alpha_3', '\alpha_4', '\alpha_5', ...
        '\alpha_6', '\alpha_7', '\alpha_8', '\alpha_9', '\alpha_{10}', '\alpha_{11}', '\alpha_{12}'})
    % xtickangles(45)
end


grid on
box on
set(gca,'linewidth',1,'Fontsize',14)
set(gca,'yscale','log')
xlabel('Parameters','Fontsize',16)
ylabel('Values in appropriate units','Fontsize',16)
if plt_save ==1
    if plt_rates ==1
    xlim([0.5 3.5])
    saveas(gcf,'D=0.01_rates.fig')
    elseif plt_all ==1
    saveas(gcf,'D=0.01_all_params.fig')  
    xlim([0.5 15.5])
    end
end

end

%% call error_all_cell_jan18 function with an internal change :  include 'flag' as 
%%%% final argument to Error_Cell_Jan18_single with D =0.01, l= 11,
%%%% rates and other parameters fixed to be mean of all iteration
%%%% corresponding to D =0.01.
sp_vec = [13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5]; 

if plt_kym ==1 
   iD = 4;
   data = load(['Optim_params_D=',num2str(D),'_tol=eminus3.dat']); 
   par_vec = data(1,:);
   %%%%% take mean 
   for icell=1:12
       dir = ['Cell',num2str(icell)]; 
       tes = Error_Cell_Jun18_single(dir,icell,sp_vec(icell),l,Dvec(iD),par_vec,dxvec(iD),'flag') 
       if plt_save ==1
       saveas(gca,['Kymograph_Cell',num2str(icell),'.fig'])
       saveas(gca,['Kymograph_Cell',num2str(icell),'.png'])
       end
   end
end 







