sp_vec = [13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5];
Dvec = [0 0.01 0.1 1];
tolf = 1e-3;  tolx= 1e-3;
dtvec = [0.1 0.1 0.05 0.01]; 
dxvec = [0.1 0.05 0.1 0.5]; 
Tot_iter = 20; 

%%%%  number of iterations with different initial conditions. 
%%%%% initial parameter values and lower-upper bounds. 
%%%%% x0 = [l koffout koffin koin alpha], these values are taken by
%%%%% suitably comparing with the first manual minimisation codes. 
iD = 4;
% for iD = 1:size(Dvec,2)
D = Dvec(iD);     
alpha0 = ones(1,12); 
x0 = [11 0.01 0.05 1 10*alpha0]; 
LB = [8 0 0 0 0*alpha0];
UB = [15 5 5 1 50*alpha0];

%%%% perform several iterations : at the end of each iteration, add noise
%%%% and start the next. 
param_iter = zeros(Tot_iter,16);
rang = [10 2 2 1 10*ones(1,12)]; 

tic 
for niter = 1:Tot_iter
    [fin_params,fval,exflag,history] = all_cell_optim(D,dxvec(iD),x0,LB,UB,tolf,tolx); 
    param_iter(niter,:)=fin_params; 
    dlmwrite(['history_D=',num2str(D),'_niter=',num2str(niter),'.dat'],history,'delimiter','\t')
    x0 = fin_params + rang.*rand(1,16); 
end
runtime = toc; 
dlmwrite(['Optim_params_D=',num2str(D),'_tol=eminus3.dat'],param_iter,'delimiter','\t')

% end

function [x,fval,exflag,history] = all_cell_optim(D,dx,x0,LB,UB,tolf,tolx)
history = [] ;
%     options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolFun',tolf...
%                             ,'TolX',tolx,'OutputFcn',@myoutput);
options = optimset('Display','iter','TolFun',tolf...
                            ,'TolX',tolx,'OutputFcn',@myoutput);
% options.MaxIter=16000;
% options.MaxFunEvals = 1000000;
[x,fval,exflag]=fminsearchbnd(@(x0) Error_all_cell_Jan18(D,dx,x0),x0,LB,UB,options);

%%%%% nested function to record history of iteration
function stop = myoutput(x,optimvalues,state)
stop = false;
  if isequal(state,'iter')
            history = [history; x];
  end 
end 

end



%%%%%  to calculate errors, use the optimised parameter values and call the
%%%%%  error function. 


% alpha0 = 0.1*ones(1,12); 
% l = 10.83;
% 
% 
% x0 = [0.02 0.01 0.05 0.1 alpha0];    
% LB = [0 0 0 0 ones(1,12)];
% UB = [5 2 2 2 2*ones(1,12)];
% 
% options = optimset('Display','iter','FunValCheck','on','MaxIter',16000,...
%     'MaxFunEvals',1000000,'PlotFcns',@optimplotfval,'TolFun',1e-4,'TolX',1e-4);
% 
% [opt_par,fval,exitflag,output]=fminsearchbnd(@(x0)...
%                                Error_all_cell_Oct17(l,x0),x0,LB,UB,options);

            
 

% %% Use opt_par to produce kymographs
% sp_vec = [13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5]; 
% err = zeros(1,12);
% for icell =1:12
%     dir = ['Cell',num2str(icell)]; 
%     sp_len = sp_vec(icell); 
%     par2 = opt_par(1:4)
%     alpha = opt_par(icell+4) 
%     
%     err(icell) =  Error_Cell_Oct17(dir,icell,sp_len,l,alpha,par2,'flag');
%     saveas(gcf,['Kymograph_Cell',num2str(icell) ,'.pdf'])
%     saveas(gcf,['Kymograph_Cell',num2str(icell) ,'.fig'])
% 
% end
% 
% f = norm(err);


% %% Sensitivity analysis 
% al_vec = opt_par(5:end); 
% 
% Daxis = 0.001:0.001:0.014;
% koff1_axis = 0.001:0.001:0.014;
% koff2_axis = 0.064:0.001:0.073;
% kon2_axis = 0.155:0.001:0.165;
% 
% for iD = 1:size(Daxis,2)
%     for ikoff1 = 1:size(koff1_axis,2)
%         for ikoff2 = 1:size(koff2_axis,2)
%             for ikon2 = 1:size(kon2_axis,2)
%              
%              par = [Daxis(iD) koff1_axis(ikoff1) koff2_axis(ikoff2) kon2_axis(ikon2)]    
%              objfun_map (iD,ikoff1,ikoff2,ikon2) = Error_all_cell_Oct17(l,[par al_vec]);
%                
%             end
%         end
%     end
% end
    

%% Second derivative 

% minim=round(opt_par(1:4),3); 
% di = 0.001;
% 
% [v1,i1] = ismember(minim(1),Daxis); 
% [v2,i2] = ismember(minim(2),koff1_axis); 
% [v3,i3] = ismember(minim(3),koff2_axis); 
% [v4,i4] = ismember(minim(4),kon2_axis); 
% 
% d2_ax1 = (objfun_map(i1+1,i2,i3,i4) + objfun_map(i1-1,i2,i3,i4)...
%             -2*objfun_map(i1,i2,i3,i4))/di/di;
% 
% 
% d2_ax2 = (objfun_map(i1,i2+1,i3,i4) + objfun_map(i1,i2-1,i3,i4)...
%             -2*objfun_map(i1,i2,i3,i4))/di/di;
% 
% d2_ax3 = (objfun_map(i1,i2,i3+1,i4) + objfun_map(i1,i2,i3-1,i4)...
%             -2*objfun_map(i1,i2,i3,i4))/di/di;
% 
% d2_ax4 = (objfun_map(i1,i2,i3,i4+1) + objfun_map(i1,i2,i3,i4-1)...
%             -2*objfun_map(i1,i2,i3,i4))/di/di;
% 
% 
% sens_vec = sqrt(objfun_map(i1,i2,i3,i4)./[d2_ax1 d2_ax2 d2_ax3 d2_ax4])









