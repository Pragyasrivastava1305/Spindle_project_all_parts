lvec=2:30; 

for il=1:size(lvec,2)
    param=lvec(il)
    errf(il) = Error_steady_state(param);
end



%%%%%%%%%%%%%% Error function %%%%%%%%%%%%%%%%
function fun = Error_steady_state(param)


% DL = params(1);
l = param;
% kon1 = params(3); kon2 = params(4);
% koff1= params(5); koff2= params(6);

kon1=60; kon2=15; koff=60; DL=1; L=36.0088; 
%  l=5; 

x=-L:0.2369:L;
I1=find(x<=-l);
I2=find(x>-l & x <=l);
I3=find(x>l);

% vel_vec=[-3 -2.5 -2 -1 -0.5 0 0.5 1 1.5 2 2.5];
 vel_vec= [-1 -0.5 0 0.5 1];


for ivel=1:size(vel_vec,2)
V=vel_vec(ivel);

lp=2*DL/(sqrt(V^2+4*DL*koff)+V); 
lm=2*DL/(sqrt(V^2+4*DL*koff)-V); 

A1=-(kon1-kon2)*lm*sinh(l/lm)*(coth(L/lm)+1)/(koff*(lp+lm));
B1=-(kon1-kon2)*lp*sinh(l/lp)*(coth(L/lp)-1)/(koff*(lp+lm));

A2=-(kon1-kon2)*lm*sinh((l-L)/lm)*csch(L/lm)/(koff*(lp+lm));
B2=-(kon1-kon2)*lp*sinh((l-L)/lp)*csch(L/lp)/(koff*(lp+lm));

A3=-(kon1-kon2)*lm*sinh(l/lm)*(coth(L/lm)-1)/(koff*(lp+lm));
B3=-(kon1-kon2)*lp*sinh(l/lp)*(coth(L/lp)+1)/(koff*(lp+lm));


%%%%  profile from theory
lgn_th(I1)=A1*exp(x(I1)/lm)+B1*exp(-x(I1)/lp)+kon1/koff;
    
lgn_th(I2)=A2*exp(x(I2)/lm)+B2*exp(-x(I2)/lp)+kon2/koff;

lgn_th(I3)=A3*exp(x(I3)/lm)+B3*exp(-x(I3)/lp)+kon1/koff;

Data=load(['mean_curve_V=',num2str(vel_vec(ivel)),'.dat']);

err1 = norm(lgn_th(I1)-Data(I1,2));
err2 = norm(lgn_th(I2)-Data(I2,2));
err3 = norm(lgn_th(I3)-Data(I3,2));


Err(ivel) = err1+err2+err3; 




 figure(1)
 hold off
  subplot(3,4,ivel)
%  subplot(2,3,ivel)
  plot(Data(:,1),Data(:,2)/Data(1,2),'-.r','linewidth',2,'Markersize',20);
  hold on
  
  plot(x,lgn_th,'k','linewidth',2,'Markersize',20)
  hold on
  drawnow

  dsize=60;
  s2=scatter(0,0); hold on; s2.Marker='o';
  set(s2,'SizeData',dsize,'MarkerFaceColor',[0 0.75 0],'MarkerEdgeColor',[0 0.75 0])
  
  title(['V=',num2str(vel_vec(ivel))],'Fontsize',16)
%   ylim([0 7000]); 
  xlim([-L L])
  drawnow
  
    
end

% close(fig)

fun = norm(Err); 
end




