% l = 10.83; 
% Dvec = [0 0.001 0.01 0.1];
ncount =0; 
Dvec = [0 0.01 0.1 1];
param_small = []; 
Tot_iter = 20; 
plt_err =0; 
plt_min =1;
plt_kym =0;
dxvec = [0.1 0.05 0.1 0.5]; 

for iD =1:size(Dvec,2)
    D = Dvec(iD);
    data = load(['Optim_params_D=',num2str(D),'_tol=eminus3.dat']); 
    
    for iter=1:Tot_iter
        [iD iter]
        param= data(iter,:);
        err(iD,iter) = Error_all_cell_Jan18(D,dxvec(iD),param);
        
    end
    if plt_err ==1
    plot(1:Tot_iter,err(iD,:),'o'); hold on ; drawnow
    end
end

%%
rc = linspace(0,1,20); 
if plt_min ==1
%%%%%  minimum appears for D=0.01, iter=7;
%%%%% generate final kymographs 
D = 0.01;
data = load(['Optim_params_D=',num2str(D),'_tol=eminus3.dat']); 
for iter =1:20
if iter ~= 7
plot(1:4,data(iter,1:4),'o','color',[rc(iter) 0 0],'Markersize',8,'linewidth',2); drawnow;  hold on   
end

% plot(1:4,data(7,1:4),'o','Markersize',8,'linewidth',2); drawnow;  hold on
end
axis square
xticks(1:4)
% xticklabels({'l', 'k^{off}_{out}', 'k^{off}_{in}','k^{on}_out', '\alpha_1' ,  '\alpha_2', '\alpha_3', '\alpha_4', '\alpha_5', '\alpha_6', '\alpha_7', '\alpha_8', '\alpha_9', '\alpha_{10}', '\alpha_{11}', '\alpha_{12}'})
% xtickangles(45)
xticklabels({'l', 'k^{off}_{out}', 'k^{off}_{in}','k^{on}_{out}'})
grid on
box on
set(gca,'linewidth',1)
set(gca,'yscale','log')
xlabel('Parameters')
ylabel('Values in appropriate units')
xlim([0.5 4.5])
end

%% call error_all_cell_jan18 function with an internal change :  include 'flag' as 
%%%% final argument to Error_Cell_Jan18_single with D =0.01, l= 15,
%%%% koff_out= 0.06, koff_in=0.19, kon_in =1
sp_vec = [13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5]; 

if plt_kym ==1 
   iD = 2;
   for icell=1:12
       dir = ['Cell',num2str(icell)]; 
       tes= Error_Cell_Jan18_single(dir,dxvec(iD),icell,sp_vec(icell),Dvec(iD),data(7,1:end),'flag') 
       saveas(gca,['Kymograph_Cell',num2str(icell),'.fig'])
       saveas(gca,['Kymograph_Cell',num2str(icell),'.png'])
   end
end 







