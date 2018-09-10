ecc= 0.2;
Phase_Matrices; 
%%%% Boundaries calculated From mathematics notebook 
if ecc == 0.2
    Math_data = load('Ph_bds_ecc02.csv'); 
else 
    Math_data = load('Ph_bds_ecc087.csv'); 
  
end
lcvec = Math_data(:,1); 
bd_min_0 = Math_data(:,2);
bd_stab_0 = Math_data(:,3);
bd_min_pihalf = Math_data(:,4);
bd_stab_pihalf= Math_data(:,5);

%%%% Run Phase_boundaries.m before running this code, as it uses the
%%%% outputs generated from Phase_boundaries.m
figure;
Lm_vec = 0.1:0.1:2;
ls_vec = 0.1:0.1:1;
  
dsize = 30; 

for ils =1:size(ls_vec,2)
both_st = find(both_yes(:,ils));
none_reach = find(both_not(:,ils)); 
st0_pi2_re =  find(z_st_pi2_re(:,ils)); 
pi2st_0_re = find(pi2_st_0_re(:,ils));
only_0reach = find(only_z_reach(:,ils)); 
only_pi2_re = find(only_pi2_reach(:,ils));
only_ob = find(only_oblique(:,ils));
ob_n_0 = find(ob_plus_0(:,ils));

%if size(vec,1)~=0
     %%%% both stable/ reachable : orange
     s1= scatter(ls_vec(ils)*ones(size(both_st,1),1), Lm_vec(both_st)); 
     hold on; s1.Marker='o';  %%% wrt lab
     set(s1,'SizeData',dsize,'MarkerFaceColor',[1 0.44 0.37],'MarkerEdgeColor',[1 0.44 0.37])
     
     hold on
     %%%% only 0 stable : pale yellow, pi/2 reachable
     s2= scatter(ls_vec(ils)*ones(size(st0_pi2_re,1),1), Lm_vec(st0_pi2_re)); 
     hold on; s2.Marker='o';  %%% wrt lab
     set(s2,'SizeData',dsize,'MarkerFaceColor',[0.94 0.86 0.51],'MarkerEdgeColor',[0.94 0.86 0.51])
     
     hold on
     %%%% only pi/2 stable, 0 reachable %%% apricot (shade of red)
     s3= scatter(ls_vec(ils)*ones(size(pi2st_0_re,1),1), Lm_vec(pi2st_0_re)); 
     hold on; s3.Marker='o';  %%% wrt lab
     set(s3,'SizeData',dsize,'MarkerFaceColor',[0.98 0.81 0.69],'MarkerEdgeColor',[0.98 0.81 0.69])
     
     hold on
     %%%% none reachable : a lighter shade of blue  
     s4= scatter(ls_vec(ils)*ones(size(none_reach,1),1), Lm_vec(none_reach)); 
     hold on; s4.Marker='o';  %%% wrt lab
     set(s4,'SizeData',dsize,'MarkerFaceColor',[0.86 0.91 0.96],'MarkerEdgeColor',[0.86 0.91 0.96])
     
     hold on
     %%%% 0 reachable  & stable 
     s5= scatter(ls_vec(ils)*ones(size(only_0reach,1),1), Lm_vec(only_0reach)); 
     hold on; s4.Marker='o';  %%% wrt lab
     set(s5,'SizeData',dsize,'MarkerFaceColor','m','MarkerEdgeColor','m')
     
     hold on
     %%%% pi/2 reachable & stable :  another shade of blue/green
     s6= scatter(ls_vec(ils)*ones(size(only_pi2_re,1),1), Lm_vec(only_pi2_re)); 
     hold on; s4.Marker='o';  %%% wrt lab
     set(s6,'SizeData',dsize,'MarkerFaceColor',[0.3 0.72 0.65 ],'MarkerEdgeColor',[0.3 0.72 0.65 ])
     
     hold on
     %%%% only oblique stable %%% dark purple
     s7= scatter(ls_vec(ils)*ones(size(only_ob,1),1), Lm_vec(only_ob)); 
     hold on; s4.Marker='o';  %%% wrt lab
     set(s6,'SizeData',dsize,'MarkerFaceColor',[0.6 0 0.6 ],'MarkerEdgeColor',[0.6 0 0.6 ])
     
     hold on
     %%%% oblqiue and 0 stable %% gray
     s8= scatter(ls_vec(ils)*ones(size(ob_n_0,1),1), Lm_vec(ob_n_0)); 
     hold on; s4.Marker='o';  %%% wrt lab
     set(s6,'SizeData',dsize,'MarkerFaceColor',[0.5 0.5 0.5 ],'MarkerEdgeColor',[0.5 0.5 0.5 ])
     
%end

xlim([0 1]);  ylim([0 2.2])

end

hold on 
%%%% plot boundaries calculated from mathematica
plot(lcvec,bd_min_pihalf,'-.','color',[0.04 0.07 0.57],'linewidth',2); 
plot(lcvec,bd_min_0,'color',[0.04 0.07 0.57],'linewidth',2); 
plot(lcvec,bd_stab_pihalf,'-.','color',[0.71 0.05 0.15],'linewidth',2); 
plot(lcvec,bd_stab_0,'color',[0.71 0.05 0.15],'linewidth',2); 


set(gca,'Fontsize',15);
box on 
grid on; axis square
xlim([0.05 0.95])
ylim([0 2.1])
xlabel('l_c/L','FontSize',16);
ylabel('l_m/L','FontSize',16); 






