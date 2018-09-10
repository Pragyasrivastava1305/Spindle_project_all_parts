function f = Error_all_cell_Jan18(D,dx,params)

%%%%%%% INPUT 
%%%% params = [l koff_far koff_near kon_near alpha_vec];
%%%% f = error

sp_vec = [13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5]; 
err = zeros(1,12);
for icell =1:12
dir = ['Cell',num2str(icell)]; 
sp_len = sp_vec(icell); 
par2 = params(1:4);
alpha = params(icell+4);
    
   if D ==0
       err(icell) =  Error_Cell_Jan18_D0(dir,dx,icell,sp_len,[par2 alpha]);
   else
       err(icell) =  Error_Cell_Jan18_single(dir,dx,icell,sp_len,D,[par2 alpha],'flag');
%        drawnow
   end
end
f = norm(err);
% if isnan(f)
% error('Nan encountered')
% end

end 
