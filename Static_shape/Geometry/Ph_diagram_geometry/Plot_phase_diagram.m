%%%%  CODE TO PLOT DIFFERENT REGIONS OF THE PHASE DIAGRAM
%%%%  AND PHASE BOUNDARIES AS CALCULATED FROM mathematica
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
ls_vec = 0.1:0.1:0.9;
  
dsize = 30; 

for ils =1:size(ls_vec,2)
vec_both = find(both_yes(:,ils));
vec_none = find(both_not(:,ils)); 
vec_0st =  find(z_stable(:,ils)); 
vec_pi2st = find(pi2_stable(:,ils));
vec_0rech = find(z_mat(:,ils)); 
vec_pi2rech = find(pi2_mat(:,ils)); 

%if size(vec,1)~=0
     %%%% both stable/ reachable : orange
     s1= scatter(ls_vec(ils)*ones(size(vec_both,1),1), Lm_vec(vec_both)); 
     hold on; s1.Marker='o';  %%% wrt lab
     set(s1,'SizeData',dsize,'MarkerFaceColor',[1 0.44 0.37],'MarkerEdgeColor',[1 0.44 0.37])
     
     hold on
     %%%% only 0 stable : pale yellow
     s2= scatter(ls_vec(ils)*ones(size(vec_0st,1),1), Lm_vec(vec_0st)); 
     hold on; s2.Marker='o';  %%% wrt lab
     set(s2,'SizeData',dsize,'MarkerFaceColor',[0.94 0.86 0.51],'MarkerEdgeColor',[0.94 0.86 0.51])
     
     hold on
     %%%% only pi/2 stable %%% apricot (shade of red)
     s3= scatter(ls_vec(ils)*ones(size(vec_pi2st,1),1), Lm_vec(vec_pi2st)); 
     hold on; s3.Marker='o';  %%% wrt lab
     set(s3,'SizeData',dsize,'MarkerFaceColor',[0.98 0.81 0.69],'MarkerEdgeColor',[0.98 0.81 0.69])
     
     hold on
     %%%% none reachable : a lighter shade of blue  
     s4= scatter(ls_vec(ils)*ones(size(vec_none,1),1), Lm_vec(vec_none)); 
     hold on; s4.Marker='o';  %%% wrt lab
     set(s4,'SizeData',dsize,'MarkerFaceColor',[0.86 0.91 0.96],'MarkerEdgeColor',[0.86 0.91 0.96])
     
     hold on
     %%%% 0 reachable     
     s5= scatter(ls_vec(ils)*ones(size(vec_0rech,1),1), Lm_vec(vec_0rech)); 
     hold on; s4.Marker='o';  %%% wrt lab
     set(s5,'SizeData',dsize,'MarkerFaceColor','m','MarkerEdgeColor','m')
     
     hold on
     %%%% pi/2 reachable :  another shade of blue/green
     s6= scatter(ls_vec(ils)*ones(size(vec_pi2rech,1),1), Lm_vec(vec_pi2rech)); 
     hold on; s4.Marker='o';  %%% wrt lab
     set(s6,'SizeData',dsize,'MarkerFaceColor',[0.3 0.72 0.65 ],'MarkerEdgeColor',[0.3 0.72 0.65 ])
     
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






