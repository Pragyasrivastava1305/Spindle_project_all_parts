Data = xlsread('Hypo_test.xlsx'); 
Cell_id = Data(:,4); 
sp_Angle_30 = 90-Data(:,1);
sp_Angle_final = 90-Data(:,2); 
sp_Angle_init = 90-Data(:,3);
Rounding_time = Data(:,6); 
rsq = Data(:,7); 
cell_wid = Data(:,8);
NEB_AR = Data(:,9); 

%%%%  find those cells for which an exponential fit to rounding dynamics 
%%%%  can be made more reliably. rsquare >  0.8. 
good_cells = find(rsq > 0.8); 
Ang_I = sp_Angle_init(good_cells);
Ang_F = sp_Angle_final(good_cells);
tau = Rounding_time(good_cells); 
wid_I = cell_wid(good_cells); 
AR_I = NEB_AR(good_cells);

%%%%% divide data in 4 groups : 
%%%% grp1 = good alignment, no substantial correction
%%%% grp2 = bad alignment, wrong correction
%%%% grp3 = bad alignment, no substantial correction
%%%% grp4 = good alignment, substantial correction
grp1 = find(Ang_I < 45 & Ang_F <  45 ); 
grp2 = find(Ang_I < 45 & Ang_F > 45 ); 
grp3 = find(Ang_I > 45 & Ang_F >  45 ); 
grp4 = find(Ang_I > 45 & Ang_F <  45 ); 


good_align = [grp1'  grp4'];
bad_align = [grp2'  grp3'];

plot(Ang_I(grp1),Ang_F(grp1),'o','linewidth',2,'Markersize',8); hold on ; drawnow
plot(Ang_I(grp2),Ang_F(grp2),'o','linewidth',2,'Markersize',8)
plot(Ang_I(grp3),Ang_F(grp3),'o','linewidth',2,'Markersize',8)
plot(Ang_I(grp4),Ang_F(grp4),'o','linewidth',2,'Markersize',8)
box on 
grid on
xlabel('Initial Angle', 'FontSize',16 )
ylabel('Final Angle', 'FontSize',16 )
set(gca, 'FontSize',14 )


%%%% group 2 and group 3 have faster rounding. check the difference in
%%% mean and intra vs inter group std. 
%%%% Rounding times. 
mTau_grp1 = mean(tau(grp1)); sTau_grp1 = std(tau(grp1)); 
mTau_grp2 = mean(tau(grp2)); sTau_grp2 = std(tau(grp2)); 
mTau_grp3 = mean(tau(grp3)); sTau_grp3 = std(tau(grp3)); 
mTau_grp4 = mean(tau(grp4)); sTau_grp4 = std(tau(grp4)); 
mTau_all = mean(tau); sTau_all = std(tau); 
mTau_vec = [mTau_grp1 mTau_grp2 mTau_grp3 mTau_grp4];
sTau_vec = [sTau_grp1 sTau_grp2 sTau_grp3 sTau_grp4];

%%%% Initial widths
mWid_grp1 = mean(wid_I(grp1)); sWid_grp1 = std(wid_I(grp1)); 
mWid_grp2 = mean(wid_I(grp2)); sWid_grp2 = std(wid_I(grp2)); 
mWid_grp3 = mean(wid_I(grp3)); sWid_grp3 = std(wid_I(grp3)); 
mWid_grp4 = mean(wid_I(grp4)); sWid_grp4 = std(wid_I(grp4)); 
mWid_all = mean(wid_I); sWid_all = std(wid_I); 
mWid_vec = [mWid_grp1 mWid_grp2 mWid_grp3 mWid_grp4];
sWid_vec = [sWid_grp1 sWid_grp2 sWid_grp3 sWid_grp4];


%%%% Initial ARs
mAR_grp1 = mean(AR_I(grp1)); sAR_grp1 = std(AR_I(grp1)); 
mAR_grp2 = mean(AR_I(grp2)); sAR_grp2 = std(AR_I(grp2)); 
mAR_grp3 = mean(AR_I(grp3)); sAR_grp3 = std(AR_I(grp3)); 
mAR_grp4 = mean(AR_I(grp4)); sAR_grp4 = std(AR_I(grp4));
mAR_all = mean(AR_I); sAR_all = std(AR_I); 
mAR_vec = [mAR_grp1 mAR_grp2 mAR_grp3 mAR_grp4];
sAR_vec = [sAR_grp1 sAR_grp2 sAR_grp3 sAR_grp4];












