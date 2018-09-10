function [th1,th2,rate] = Rates_dynamic_circular_DNA_May18(x,th,s,a,b,dnaR,kn1,kn2,varargin)
%%%%% INPUT variables   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% x = grid vector on which LGN concentration is defined. has different
%%%%% spacing from s which has uniform spacing in theta but not in length. 
%%%%% Ranges of the two are same. th = theta vector, uniform spacing. 
%%%%% (a,b) = semi-major and minor axes of ellipse, ang =  angle of DNA 
%%%%% disc wrt long axis; lA = length of DNA arm,lR = Ran disc radius,
%%%%% (kn1,kn2) = rates outside, and inside inhibited region.
%%%%% OUTPUT variables: th1,2 = the two inhibited regions, can be an empty
%%%%% vector when there is no crossing, or vectors of length 2 or 4. 
%%%%% rate = on/off rates. To help with troubleshooting: uncomment plot
%%%%% commands. 

%%%% DNa is circular in this case and the only difference from dna
%%%% non-circular is that the distances are measured from the origin. 

%   ang = round(-1.55+pi/2,3);  
%   lA = ArmR; lR=RanR; kn1 = 0.5; kn2=1; 
xl_el =-40; yl_el = -40; 
xh_el =40; yh_el = 40; 

vl = 0.8;
flag = ismember('flag',varargin); 
lw = 2; ms =3;
r=@(th)a./sqrt(a^2/b^2 - (a^2/b^2-1)*cos(th).^2);
xunit = r(th).* cos(th) ;
yunit = r(th).* sin(th) ;

% % %%%% compatibility checks 
%  if dnaR >= a 
%      error('DNA is too big: bigger than cell length')
%  end

dist_PtoN1 = sqrt((r(th).*cos(th)).^2+(r(th).*sin(th)).^2);
dist_PtoN2 = sqrt((r(th).*cos(th)).^2+(r(th).*sin(th)).^2);
%%%%  location of the boundaries between regions of kn1 and kn2.
Ivec1 = find(abs(diff(dnaR - dist_PtoN1 > 0))); 
Ivec2 = find(abs(diff(dnaR - dist_PtoN2 > 0))); 

%%%% When the circles exactly touch the ellipse, sometimes 
%%%% size(Ivec1,2)= and size(Ivec2,2) =2. To rectify this we add following
%%%% lines of code. Same error at one grid point may happen when
%%%% size(Ivec1,2) ==2 while Ivec2=[] or vice versa. 

if (size(Ivec1,2)>2 && Ivec1(4) == size(th,2)-1 && Ivec1(1) ==1)
    Ivec1 = [Ivec1(2) Ivec1(3)];
end 

%%% correct for only one grid point intersection. 
if (size(Ivec1,2)==2 && abs(diff(Ivec1))==1)
    Ivec1 = [];
end 

%%% correct for only one grid point intersection. 
if (size(Ivec2,2)==2 && abs(diff(Ivec2))==1)
    Ivec2 = [];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rate = kn1*ones(1,size(x,2)); 


%%%% 0 CROSSING CASE
if (size(Ivec1,2)<1)     
%   warning('Ran range does not reach cell periphery: no inhibition')    
    if dnaR < a
    th1 = [];
    th2 = [];
    elseif dnaR >= a
    th1 = [];
    th2 = [];  
%     th1 = [th(1) th(end)]; 
%     th2 = [th(1) th(end)]; 
    rate = kn2*ones(1,size(x,2));     
    end
    
    %%%% plot if flag is passed 
    if flag==true
    if dnaR < a    
    plot(xunit,yunit,'ro','Markersize',ms);
%     hold on; plot([N1x,N2x],[N1y,N2y],'og-','linewidth',lw,'Markersize',3*ms...
%           ,'Markeredgecolor','k','Markerfacecolor','g'); 
    elseif dnaR >= a
    plot(xunit,yunit,'go','Markersize',ms);    
    end
    
    hold on; 
    circle(0,0,dnaR,lw,'g');   %%% ran discs 
    circle(0,0,dnaR,lw,'g');  
    xlim([xl_el  xh_el]); ylim([yl_el,yh_el])
    axis equal tight
    hold off
    end
    
%%%% 2 CROSSINGS CASE
elseif (size(Ivec1,2)==2 )
  
    %%%% inhibited region 1
    s1 = [s(Ivec1(1)+1) s(Ivec1(end))] ;
    th1 = [th(Ivec1(1)+1) th(Ivec1(end))] ;
    
    %%%%  set up a condition to choose correct interval 
    if (dnaR-dist_PtoN1(Ivec1(1)+2)) > 0  
        Exch1 =0;
    else
        Exch1 =1;
    end
    
    
    if Exch1 == 0
       rate(x >= s1(1) & x <= s1(2)) = kn2; 
       
    elseif (Exch1 ==1)
        th1 = [th1(2) th1(1)]; %%%% interchange  so the region is between th1(2) and th1(1)
        s1 = [s1(2) s1(1)];
%         rate = kn1*ones(1,size(x,2)); 
        rate(x >= s1(1)) = kn2; 
        rate(x <= s1(2)) = kn2; 
        %plot(rate)
    end
    
    %%%%  set up a condition to choose correct interval 
    if (dnaR-dist_PtoN2(Ivec2(1)+2)) > 0  
        Exch2 =0;
    else
        Exch2 =1;
    end
    
    s2 = [s(Ivec2(1)+1) s(Ivec2(end))] ;
    th2 = [th(Ivec2(1)+1) th(Ivec2(end))] ;
    
    if Exch2 == 0
       rate(x >= s1(1) & x <= s1(2)) = kn2; 
    elseif Exch2 == 1
       th2 = [th2(2) th2(1)];
       s2 = [s2(2) s2(1)]; 
       rate(x >= s2(1)) = kn2; 
       rate(x <= s2(2)) = kn2; 
    end
    
     
    reg1 = (th >= th1(1) & th < th1(2));  %%% inhibited region 1
    reg2 = (th >= th2(1) & th < th2(2));  %%% inhibited region 2

    if (th1(1) > th1(2))
        reg1 = ( th <= th1(2) | th >= th1(1));
    end
    
    if (th2(1) > th2(2))
        reg2 = ( th <= th2(2) | th >= th2(1));
    end
    
    %%%%% plot if flag is passed
    if flag == true
    plot(xunit,yunit,'ro','Markersize',3); hold on; 
    plot(xunit(reg1),yunit(reg1),'go',xunit(reg2),yunit(reg2),'go','Markersize',3);
    hold on; circle(0,0,dnaR,lw,'g'); circle(0,0,dnaR,lw,'g');
    xlim([-20 20]); ylim([-20,20])

    axis equal tight
    hold off   
    end
    
%%%% 4 CROSSINGS CASE
elseif (size(Ivec1,2)>2 )
    s1 = [s(Ivec1(1)+1) s(Ivec1(2))  s(Ivec1(3)+1)  s(Ivec1(4))] ;
    s2 = [s(Ivec2(1)+1) s(Ivec2(2))  s(Ivec2(3)+1)  s(Ivec2(4))] ;
    
    th1 = [th(Ivec1(1)+1) th(Ivec1(2))  th(Ivec1(3)+1)  th(Ivec1(4))] ;
    th2 = [th(Ivec2(1)+1) th(Ivec2(2))  th(Ivec2(3)+1)  th(Ivec2(4))] ;
    
    rate(x >= s2(1) & x <= s2(2)) = kn2; 
    rate(x >= s2(3) & x <= s2(4)) = kn2; 
    
    reg3 = (th >= th2(1) & th < th2(2));  %%% inhib
    reg4 = (th >= th2(3) & th < th2(4));  %%% inhib
    
    %%%% since part 1 or part 2 of the region 1 can cross the linear
    %%%% boundaries, the order of angles should be changed to pick up
    %%%% correct intervals. Define two conditions : 
    %%% condition that region 1 crosses left linear boundary and is in
    %%% first and 4 angular quad : only one of the points will be in upper
    %%% half 
    
    cond1 = (th1(1) >=0 && th1(1) < pi) && ...
           (th1(2) >=pi && th1(2)<= 2*pi) && ... 
           (th1(3) >=pi && th1(3)<= 2*pi) && ...
           (th1(4) >=pi && th1(4)<= 2*pi); 
   
    %%% Same as above for part 2 of region 1
    
    cond2 = (th1(1) >=0 && th1(1) < pi) && ...
           (th1(2) >=0 && th1(2)<= pi) && ... 
           (th1(3) >=0 && th1(3)<= pi) && ...
           (th1(4) >=pi && th1(4)<=2*pi); 

   %%% correct for the order of angles
     if cond1 
         th1 = [th1(4) th1(1) th1(2) th1(3)];
         s1 = [s1(4) s1(1) s1(2) s1(3)]; 
         
         rate(x >= s1(1) | x <= s1(2)) = kn2; 
         rate(x >= s1(3) & x <= s1(4)) = kn2; 
         
         reg1 = (th >= th1(1) | th < th1(2));  %%% inhib
         reg2 = (th >= th1(3) & th < th1(4));  %%% inhib
         
     elseif cond2
         th1 = [th1(2) th1(3) th1(4) th1(1)];
         s1 = [s1(2) s1(3) s1(4) s1(1)]; 
          
         rate(x >= s1(1) & x <= s1(2)) = kn2; 
         rate(x >= s1(3) | x <= s1(4)) = kn2; 
          
         reg1 = (th >= th1(1) & th < th1(2));  %%% inhib
         reg2 = (th >= th1(3) | th < th1(4));  %%% inhib
     else
                  
         rate(x >= s1(1) & x <= s1(2)) = kn2; 
         rate(x >= s1(3) & x <= s1(4)) = kn2; 
         reg1 = (th >= th1(1) & th < th1(2));  %%% inhib
         reg2 = (th >= th1(3) & th < th1(4));  %%% inhib
         
     end
      
    
     
    %%%%% plot if flag is true
    if flag == true
    plot(xunit,yunit,'ro','Markersize',3); hold on; 
    plot(xunit(reg1),yunit(reg1),'go',xunit(reg2),yunit(reg2),'go',xunit(reg3),yunit(reg3),'go',...
         xunit(reg4),yunit(reg4),'go','Markersize',3);

    hold on; circle(0,0,dnaR,lw,'g'); circle(0,0,dnaR,lw,'g');
    xlim([-20 20]); ylim([-20,20])

    axis equal tight
    hold off
    end
    
end




