D=1; dx=0.5;  l=11;
param = [0.0336 0.1755 0.9999 8.8381 19.4850 18.2150 14.4760 14.0960 17.8800 9.9861 28.9140...
    23.0880 18.0630 15.7070 10.9170];

sp_vec=[13.5 15 13.5 16 13.5 13 12.5 12.5 13.5 14.5 14.5 14.5];

plt_save =0 ; 
plt_grp =1; 
grp1 = [1,3,4,5,7,10,11,12];
vec = linspace(0,1,size(grp1,2)); 
%%%%  run the for loop only on group 1 if plt_grp =1
% for igrp =1:size(grp1,2)
% icell = grp1(igrp);

for icell =1:12
dir = ['Cell',num2str(icell)]; 
sp_len = sp_vec(icell); 

%%%%  experimental LGN 
lgn_init = load(fullfile(dir,['Cell',num2str(icell),'_lgn_init.dat']));  
xlab = lgn_init(:,1);  
vec = diff(xlab);
dxlab =  vec(1);

%%%% this sets system size to be an integer multiple of numerical gridsize 
L=floor(xlab(end)/dx)*dx; 
x=(0:dx:L); N=size(x,2); 
    
%%%%%%% LGN in space time, smoothening and interpolation on numerical grid
lgn_data=load(fullfile(dir,['Cell',num2str(icell),'_conc.dat']));
Nstep=size(lgn_data,1); 

for ifr=1:Nstep
%         lgn_aver(ifr,:)=movmean(lgn_data(ifr,:),10);
        lgn_exp(ifr,:)=interp1(xlab,lgn_data(ifr,:),x);
end  
  
%%%% simulated LGN
[c_sim, err(icell)] = React_diff_solution(dir,icell,sp_len,l,D,param,dx); 

tvec = 0:3:3*size(c_sim,2)-1; 
mean_sim = dx*3*sum(sum(c_sim))/tvec(end)/L; 
mean_exp = dx*3*sum(sum(lgn_exp))/tvec(end)/L;

for ifr=1:size(tvec,2)
    Tot_sim(ifr) = sum(c_sim(:,ifr))*dx/mean_sim;  
    Tot_exp(ifr) = sum(lgn_exp(ifr,:))*dx/mean_exp;
end
figure;

plot(tvec,Tot_exp,'o-.r','linewidth',1.5); drawnow;  hold on;
%figure(2)
plot(tvec,Tot_sim,'k','linewidth',1.5); drawnow; hold on;
set(gca,'fontsize',16)
xlabel('Time (min)','fontsize',16)
ylabel('Normalised concentrations','fontsize',16)
legend('Experiments','Simulated')
xlim([0 tvec(end)+3])

axis square; grid on; box on;
if plt_save ==1
    saveas(gcf,['Compare_tot_LGN_Cell=',num2str(icell),'.png'])
    
end

clear Tot_sim Tot_exp c_sim c_exp lgn_data lgn_exp
end












