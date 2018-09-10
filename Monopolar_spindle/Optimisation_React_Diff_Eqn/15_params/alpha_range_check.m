D = 0.01; koff_out =0.07; koff_in =0.23; kon_in = 1; 
alpha_vec = [1 5 10 50 100 500 1000 5000 10000]; 
sp_vec = [13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5]; 
l = 11;
err_alpha = zeros(12,size(alpha_vec,2)); 
param = [koff_out koff_in kon_in]; 

for icell =1:12
    dir =['Cell',num2str(icell)];
    icell 
    for ial =1:size(alpha_vec,2)
        err_alpha(icell,ial) = Error_Cell_Jun18_single(dir,icell,sp_vec(icell),l,D,[param alpha_vec(ial)],0.02);
    end
    
    [almin(icell),indmin(icell)] = min(err_alpha(icell,:)); 
    plot(alpha_vec,err_alpha(icell,:),'o-.'); hold on; drawnow
end

Error_all_cell_Jun18(l,D,[param almin],0.05)