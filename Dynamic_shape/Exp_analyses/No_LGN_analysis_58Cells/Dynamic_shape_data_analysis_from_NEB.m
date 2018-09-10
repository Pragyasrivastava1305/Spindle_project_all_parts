% clear all;
% close all; 
%%%%%%% load the data in a separate file. 
% Load_dynamic_data; 
%%%%% data is extracted using the following code and the quantities are
%%%%% defined. 
Load_dynamic_data_All;
%%%%%%%
plt_save =0; 
indiv_fit =0;

loc = [0 find(diff(Dynamic_shape_all.cell_id))' size(Dynamic_shape_all,1)]; 
prompt = 'Enter the quantity to be plotted (fmt: Dynamic_shape_all.quantity): '; 

prompt2 = 'Enter 0 for raw plots, 1 for normalised plots: ';
quant_y = input(prompt); 
norm_opt = input(prompt2); 


%%%%  create a matrix for given number of cells and 33/18 times points. 33/18 time points are
%%%%  decided by once running the code to see which cell has been tracked
%%%%  longest.  Run the code with 0 as input to prompt2 and comment all
%%%%  lines containing Lmat and Lvec. 33/18 is maximum of tsize vector. 

Lmat = zeros(size(loc,2)-1,33);
color_vec = linspace(0,1,size(loc,2)-1);


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
    quant_I(icell) = quant_dyn(t0); 
    quant_F(icell) = quant_dyn(end); 
    
    tsize(icell) = size(time,1)-t0+1;
    
    
    if norm_opt == 0 
    figure(1)
    plot(time(t0:end),quant_dyn(t0:end),'.-.','color',[0.5 1-color_vec(icell) 0],'linewidth',1.2,'Markersize',16)
%     plot(time(t0:end),quant_dyn(t0:end),'o','color',[0.5 0.5 0],'linewidth',2,'Markersize',16)
    Lvec(1:tsize(icell)) = quant_dyn(t0:end); 
    Lmat(icell,:) = Lvec; 
       
    elseif norm_opt ==1
    %figure;
    norm_y = (quant_dyn(t0:end)-quant_dyn(end))/(quant_I(icell)-quant_F(icell)); 
    plot(time(t0:end),norm_y,'.-.','color',[0.5 1-color_vec(icell) 0],'linewidth',1.5)
    %plot(time(t0:end),norm_y,'-.or','linewidth',1.5,'Markersize',8)
    hold on
    
    %%%%% exponential fit to each curve to get distribution of
    %%%%% characteristic time scales. 
    if indiv_fit ==1
    figure;
    [indiv_q1, indiv_q2] = createFit(time(t0:end),norm_y);
    tau_ind(icell) = 1/abs(indiv_q1.b);
    rsq(icell) = indiv_q2.rsquare;
    ydata = indiv_q1.a*exp(indiv_q1.b*time(t0:end));
    plot(time(t0:end), ydata,'-k','linewidth',2)
    axis square
    set(gca,'fontsize',14)
    xlabel('time','fontsize',16);
    ylabel('Normalised Aspect Ratio','fontsize',16);
    
    if plt_save==1
        saveas(gca,['Cell=',num2str(icell),'.png'])    
    end
    hold on;
    drawnow
    close(gcf)
    end
    
    Lvec(1:tsize(icell)) = norm_y; 
    Lmat(icell,:) = Lvec; 
  
    end
    
    axis square
    set(gca,'fontsize',14)
    xlabel('time','fontsize',16);
    ylabel('Normalised Aspect Ratio','fontsize',16);
    
    
    hold on
    clear time 
    clear quant_dyn Lvec 
end

tvec = linspace(0,3*(max(tsize)-1),max(tsize));

for it = 1:size(tvec,2)
     I = find(Lmat(:,it)); 
     size_pnt(it) = size(I,1);
     mean_quant(it) = mean(Lmat(I,it));
     std_quant(it) = std(Lmat(I,it)); 
end

% errorbar(tvec,mean_quant,std_quant,'k','linewidth',2)
% set(gca,'fontsize',14)
% xlabel('time','fontsize',16);
% ylabel('Normalised cell area', 'fontsize',16)
% xlim([0 tvec(end)])


h1 = errorbar(tvec,mean_quant,std_quant,'k-.','linewidth',2)
hold on 

if norm_opt ==1
%%% fit to overall mean
[q1,q2] = createFit(tvec,mean_quant); 
tau = 1/abs(q1.b); 
hold on;
ydata = q1.a*exp(q1.b*tvec);
h2 = plot(tvec, ydata,'r','linewidth',2)

elseif norm_opt ==0
tnew = tvec(1:18);
quant_new = mean_quant(1:18)-mean_quant(18);
[q1,q2] = createFit(tnew,quant_new); 
ydata = q1.a*exp(q1.b*tnew)+mean_quant(18);
tau = 1/abs(q1.b);
h2 = plot(tnew, ydata,'r','linewidth',2)
hold on;  

% dlmwrite('DNA_width_fit.csv',[tnew' ydata'])
end

set(gca,'fontsize',14)
xlabel('time','fontsize',16);
% ylabel('Normalised cell area', 'fontsize',16)
xlim([0 tvec(end)])
%legend(gca,'exp(-b*t)','Data', 'Location', 'NorthEast' );

% 
% %%%%%  Sort cells from fastest to slowest 
% if norm_opt ==1
% [Srate,Irate] = sort(tau_ind); 
% end

legend([h1 h2],{'Mean curve',['Fit =',num2str(q1.a), '*exp(-t/',num2str(round(1/abs(q1.b))),')']})




%%%%%% equation for cell length 
%%%% cell_length(t) = 22.83 + q1.a*exp(q1.b*t);
%%%% [q1.a,q1.b] = [14.84, -0.24]

%%%%%% equation for cell width 
%%%% cell_width(t) = 22.35 + q1.a*exp(q1.b*t);
%%%% [q1.a,q1.b] = [1.62, -0.57]



























