close all
par_data = load('history_D=0.01_niter=1_dx=02.dat');
final_params = par_data(end,:); 
l =11; 
D = 0.01;
sp_vec = [13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5]; 
dx = 0.02;

for icell=1:12
    dir = ['Cell',num2str(icell)];
    sp_len = sp_vec(icell); 
    %set(gcf,'Renderer','painters')
    Error_Cell_Jun18_single(dir,icell,sp_len,l,D,final_params,dx,'flag')
    drawnow
    saveas(gcf,['Cell',num2str(icell)],'pdf')
    %close(gcf)
end