Dvec=[0 0.01 0.1 1];
iD = 4
D = Dvec(iD);
figure;
params_all_iter = zeros(16,10);
% for
    iter = 1
    Data = load(['history_D=',num2str(0),'_niter=',num2str(iter),'.dat']); 
    params_all_iter(:,iter) = Data(end,:); 
    scatter(1:16,params_all_iter(1:16,iter)); hold on; drawnow
%     err(iD,iter) = Error_all_cell_Jan18(D,params_all_iter(:,iter)'); 
    
% end
% 
% lm(iD) = mean(params_all_iter(1,:));
% koffo_m(iD) = mean(params_all_iter(2,:));
% koffi_m(iD) = mean(params_all_iter(3,:));
% kon_m(iD) = mean(params_all_iter(4,:));