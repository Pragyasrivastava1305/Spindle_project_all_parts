% clear all;
% close all; 
%%%%%%% load the data in a separate file. 
% Load_dynamic_data; 
%%%%% data is extracted using the following code and the quantities are
%%%%% defined. 
Load_dynamic_data_All;
%%%%%%%
plt_save =0; 

loc = [0 find(diff(Dynamic_shape_all.cell_id))' size(Dynamic_shape_all,1)]; 
prompt = 'Enter the quantity to be plotted (fmt: Dynamic_shape_all.quantity): '; 

prompt2 = 'Enter 0 for raw plots, 1 for normalised plots: ';
quant_y = input(prompt); 
norm_opt = input(prompt2); 


%%%%  create a matrix for given number of cells and 33/18 times points. 33/18 time points are
%%%%  decided by once running the code to see which cell has been tracked
%%%%  longest.  Run the code with 0 as input to prompt2 and comment all
%%%%  lines containing Lmat and Lvec. 33/18 is maximum of tsize vector. 

Lmat = zeros(size(loc,2)-1,38);
color_vec = linspace(0,0.25,size(loc,2)-1);


for icell = 1:size(loc,2)-1
    Lvec = zeros(1,size(Lmat,2));
    %%%% for check 
%   [Dynamic_shape.id(loc(icell)+1) Dynamic_shape.id(loc(icell+1)+1)] 
    
    time = Dynamic_shape_all.time(loc(icell)+1:loc(icell+1));
    
    %%% find the 0 of time.
    t0 = find(time==0); 
    zero(icell) = t0; 
    quant_dyn = quant_y(loc(icell)+1:loc(icell+1));
    
    %%%% initial and final values of the quantity
    quant_I(icell) = quant_dyn(1); 
    quant_F(icell) = quant_dyn(end); 
    
    tsize(icell) = size(time,1);
    
    %%%%% set norm_opt = 1 for normalised plots. 
    figure(1)
    if norm_opt == 0       
    plot(time,quant_dyn,'.-.','color',[0.5 0 1-color_vec(icell)],'linewidth',1.2,'Markersize',16)
    hold on; drawnow
    Lmat(icell,1:size(quant_dyn,1)) = quant_dyn;
    
    elseif norm_opt ==1
    norm_y = (quant_dyn(1:end)-quant_dyn(end))/(quant_I(icell)-quant_F(icell)); 
    plot(time,norm_y,'.-.','color',[0.5 1-color_vec(icell) 0],'linewidth',1.2,'Markersize',8); drawnow
    hold on
    Lmat(icell,1:size(quant_dyn,1)) = norm_y;
    
    end
    
    if plt_save==1
        saveas(gca,['Cell=',num2str(icell),'.png'])    
    end
    
    axis square
    set(gca,'fontsize',14)
    xlabel('time','fontsize',16);
    ylabel('Normalised Aspect Ratio','fontsize',16);
    
    hold on
    clear time 
    clear quant_dyn Lvec 
end



%%%%% plot mean and std curves 
fin_frame = 19; 
tvec = 3*(-5:max(tsize)-6); 
%tvec = linspace(,3*(max(tsize)-1),max(tsize));
%for it = 1:size(tvec,2)  %%%%  uncomment this for loop to plot upto the
                          %%%% last time point of the largest duration measurement.
for it =1:fin_frame
     I = find(Lmat(:,it)); 
     size_pnt(it) = size(I,1);
     mean_quant(it) = mean(Lmat(I,it));
     std_quant(it) = std(Lmat(I,it)); 
end

h1=errorbar(tvec(1:fin_frame),mean_quant,std_quant,'ro','linewidth',2,'markersize',8);
hold on 

tnew  = tvec(1:fin_frame);

%%% tanh fit to overall mean
[xData, yData] = prepareCurveData(tnew, mean_quant );
ft = fittype( 'a*(1-tanh((x-C)/tau))+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.913375856139019 0.980883206682952 0.360377959039892 0.63235924622541];
[fit_res,fit_opt ] = fit( xData, yData, ft, opts );

a = fit_res.a;
C = fit_res.C;
tau =  fit_res.tau;
b = fit_res.b;

ydata = a*(1-tanh((tnew-C)/tau))+b;
hold on;
h2 = plot(tnew, ydata,'k','linewidth',2);

xlim([tnew(1) tnew(end)]); 
box on; grid on

%legend([h1 h2],{'Mean curve',['Fit =',num2str(q1.a), '*exp(-t/',num2str(round(1/abs(q1.b))),')']})


%%%%%% equation for cell length 
%%%% cell_length(t) = 22.83 + q1.a*exp(q1.b*t);
%%%% [q1.a,q1.b] = [14.84, -0.24]

%%%%%% equation for cell width 
%%%% cell_width(t) = 22.35 + q1.a*exp(q1.b*t);
%%%% [q1.a,q1.b] = [1.62, -0.57]



























