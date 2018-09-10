function f = Error_all_cell_Jun18(l,D,params,dx)

%%%%%%% INPUT 
%%%% l = inhibition range 
%%%% params = [D koff_far koff_near kon_near alpha_vec];
%%%% f = error


sp_vec = [13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5]; 
err = zeros(1,12);

for icell =1:12

dir = ['Cell',num2str(icell)]; 
sp_len = sp_vec(icell); 
par2 = params(1:3);
alpha = params(icell+3);
    
   if D ==0
       err(icell) =  Error_Cell_Jun18_D0(dir,icell,sp_len,l,[par2 alpha]);
   else
       err(icell) =  Error_Cell_Jun18_single(dir,icell,sp_len,l,D,[par2 alpha],dx,'flag');
   end
% saveas(gcf,['Cell=',num2str(icell),'D=0.01.png'])   
end

f = norm(err);
end 








