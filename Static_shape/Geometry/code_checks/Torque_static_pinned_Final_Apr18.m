function [Omega,torq_dens1,torq_dens2,ang1,ang2] = Torque_static_pinned_Final_Apr18(a,b,th,s,dth,rn,rs,sp_ang,lsp,Lm,x,sys_L,conc,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  (a,b) = semi-major and minor axes of ellipse, th = angular grid
%%%%%  s = curvilinear coordinate: irregularly spaced grid.
%%%%%  dth = angular grid spacing.
%%%%%  rn,rs =  rounding decimal place for angle and position.
%%%%%  sp_ang = instantaneous angle of the spindle
%%%%%  lsp = half length of the spindle
%%%%%  Lm = range of MT, 
%%%%%  x = uniform space grid along periphery
%%%%%  L = system size
%%%%%  conc = concentration of LGN
%%%%%  varargin = 'flag' to generate figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Omega =  resultant torque on the spindle. 2]
%%%%%  [Torq_dens1, torq_dens2] = Torque density along s. 
%%%%%  [ang1, ang2] =  angular intersection around centrosome 1 and 2. 
%%%%%  Each of ang1 and ang2 can be vectors of lengths 0, 2 and 4 for zero
%%%%%  crossings, 2 crossings and 4 crossings. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%  ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  a) Calculate the range of overlap between MT and periphery.
%%%%%  Can have 3 cases :  no overlap, circle crosses periphery twice or 4
%%%%%  times. Depending on which quadrants the intersections fall into, the
%%%%%  order of angles in the interesection vectors needs to be changed to
%%%%%  ensure that right interval is picked up for torque integration.
%%%%%  Two parameters are assigned to keep track of 1.changing the order of
%%%%%  angles 2. Taking care of the linear crossing, denoted by Exch and
%%%%%  Lincross. The two happen at the same time so keeping track of one
%%%%%  parameter is enough. Exch = 1. 
%%%%%  b) interpolate conc from x to s grid 
%%%%%  c) Integrate this conc with appropriate prefactors, to calculate
%%%%%  torque on each centrosome. Care about the periodicity and wrapping
%%%%%  back of the concentration from linear to elliptical. 
%%%%%  To help with the troubleshoot: uncomment plot commands. 

%%%%%% BECAUSE OF THE SYMMETRY OF THE PROBLEM, AND IN ABSENCE OF NOISE, IT
%%%%%% IS ENOUGH TO CONSIDER THE ANGLE OF CENTROSOME 1 IN QUAD 1 WHICH ALSO
%%%%%% CONSTRAINS CENTROSOME 2 TO BE IN QUAD 3. IF NOISE IS TO BE INCLUDED,
%%%%%% IT WILL BECOME NECESSARY TO TAKE INTO ACCOUNT THE CENTROSOME 1 IN
%%%%%% 4TH QUAD AND CENTROSOME 2 IN 2ND QUAD. DUE TO THE SAME REASON, ONLY
%%%%%% NEGATIVE REGION OF CENTROSOME 1 AND, ONLY POSITIVE REASON OF 
%%%%%% CENTROSOME 2 CAN CROSS LINEAR BOUNDARY

% b=10;  a=20;
% rn=8; rs=1;  sp_ang = round(0.1071,rn); lsp = 7; Lm = 13;
% th = round(linspace(0,2*pi,360*20),rn);
% 
% dth = round(th(2)-th(1),3); 
% conc = ones(1,size(x,2));
lw =2; ms=3; vl =0.8; 
% xl_el = -20; xh_el = 20; yl_el = -20; yh_el =20; 
 
 
flag = ismember('flag',varargin); 
torq_dens1 = zeros(1,size(s,2));
torq_dens2 = zeros(1,size(s,2));
%%%%% compatibility checks 
if lsp >= b 
    error('spindle is bigger than or equal to the cell width')
end


%%%%% define shape of the ellipse and the vector along periphery.
r=@(th)a./sqrt(a^2/b^2 - (a^2/b^2-1)*cos(th).^2);
xunit = r(th).* cos(th) ;
yunit = r(th).* sin(th) ;

%%%%% angles of the two centrosomes.
%%%%% define the centrosome which is in 1st/4th quad to be cen1 and the one
%%%%% in 2nd/3rd quad to be cen2. 
if ( cos(sp_ang)>0 )
    phi1 = sp_ang; phi2 = mod(sp_ang + pi,2*pi); 
elseif ( cos(sp_ang)<0 )
    phi1 = mod(sp_ang + pi,2*pi); phi2 = mod(sp_ang,2*pi);
end

if ( mod(sp_ang,pi/2) <= dth && sp_ang~=0 )
    phi1 = pi/2; phi2 = 3*pi/2; 
end

%%% coordinates of centrosomes
C1x = lsp*cos(phi1); C1y = lsp*sin(phi1);
C2x = lsp*cos(phi2); C2y = lsp*sin(phi2);

fun1 = @(mu) Lm-sqrt(lsp^2 - 2*lsp*r(mu).*cos(mu-phi1)+ r(mu).^2);
fun2 = @(mu) Lm-sqrt(lsp^2 - 2*lsp*r(mu).*cos(mu-phi2)+ r(mu).^2);
Ivec1 = find(abs(diff( fun1(th) >= 0))); 
Ivec2 = find(abs(diff( fun2(th) >= 0))); 


%%%% When the circles exactly touch the ellipse, sometimes 
%%%% size(Ivec1,2)= and size(Ivec2,2) =2. To rectify this we add following
%%%% lines of code
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


%%%%% consider the 3 cases of Ivec1, and Ivec2.

if flag==true
%      figure(2)
    plot(xunit,yunit,'ro','Markersize',ms);
    hold on; plot([C1x,C2x],[C1y,C2y],'sq-','Color','[1 0.6 0]','linewidth',lw,'Markersize',4*ms...
         ,'Markeredgecolor','k','Markerfacecolor','[1 0.6 0]');  
    hold on; 
    circle(C1x,C1y,Lm,lw,[vl vl vl]);   %%% ran discs 
    circle(C2x,C2y,Lm,lw,[vl vl vl]);  
%     xlim([xl_el  xh_el]); ylim([yl_el,yh_el])
    axis equal tight
    drawnow
end


if (size(Ivec1,2)<1)   %%%%  THIS IF LOOP ENDS ONLY AT THE VERY END OF CODE.
    
    %%%%%%%%%%%%%%%%%%% ZERO CROSSING CASE %%%%%%%%%%%%%%%%%%%%%%%%
    warning('Astral MTs do not intersect the periphery. resulting torque is 0')
    Omega = 0;
    ang1 = [];
    ang2 = [];
    %%% torq_dens remains 0. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% TWO CROSSING CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%--------------------------------------------------------%%%%%%%%%%%%
%%%% Calculate torque on first centrosome whose overlap region with
%%%% periphery can cross the linear system's boundaries. Check for this
%%%% is whether ang1(1)<phi1<ang1(2) and same for phi2. 
    
elseif (size(Ivec1,2) ==2)
    
    %%%%%%%% CENTROSOME 1
    ang1 = round([th(Ivec1(1)+1) th(Ivec1(end))],rn); 
    
    %%% the following condition determines which angular interval should be
    %%% taken for integration (ang1(1) ang1(2)) or (ang1(2) ang1(1)). 
    if fun1(th(Ivec1(1)+2)) > 0  
        Exch1 =0;
    else
        Exch1 =1;
    end
    
    if Exch1 ==1
    [ang1] = [ang1(2) ang1(1)]; 
    end
    
    %%%%%%%%%%%% CASE 1 :  NO CROSSING OF BOUNDARIES %%%%%%%%%%%%%%%%%%%%%%
    if (Exch1 == 0)  
        %%%% input range for interpolation on s
        s1 = [s(Ivec1(1)+1) s(Ivec1(end))];
        
        %%%% range and indices on x vector
        xran1 = find(x >= s1(1) & x <= s1(2)); 
        Il1 = xran1(1)-1 ; Ir1 = xran1(end)+1; 
        
        %%%% specification in case xran1 end points are exactly on one of
        %%%% the edges of linear system.
        if Il1 ==0;  Il1 = 1; end
        if Ir1 == size(x,2)+1; Ir1 = size(x,2); end
        
        %%%%  angular and s overlap regions around centrosome 1 and 2
        arange1 = th(th >= ang1(1) & th <= ang1(2));
        srange1 = s(th >= ang1(1) & th <= ang1(2));
        
        %%%% interpolate concentration on s 
        c1interp = interp1(x(Il1:Ir1),conc(Il1:Ir1),srange1); 
        chk1 = find(isnan(c1interp));
        
        if (size(chk1,2)>=1)
            error('Interpolation error in 2 crossings, centrosome 1, case 1')
        end

        den1 = sqrt((r(arange1).*cos(arange1)-lsp*cos(phi1)).^2+...
                       (r(arange1).*sin(arange1)-lsp*sin(phi1)).^2);
              
        torq_dens1(th >= ang1(1) & th <= ang1(2)) = lsp*r(arange1).*c1interp.*sin(arange1-phi1)./den1; 
        torq1 = trapz(dth*lsp*r(arange1).*c1interp.*sin(arange1-phi1)./den1);
%         torq1 = trapz(dth*torq_dens(arange1)); 
        
    %%%%%%%%%%%%%%%%%%% CASE 2 : REGION CROSSES LINEAR BOUNDARIES:
    %%%%%%%%%%%%%%%%%%% Correct order of limits is already obtained. 
    %%%%%%%%% The negative region needs to be broken in two while
    %%%%%%%%% calculation from positive region remains same.         
    elseif (Exch1 == 1)        
        subr1 = [s(Ivec1(end)) x(end)];    %% region in 4th quad 
        subr2 = [0 s(Ivec1(1)+1)];         %% region 1st quad        
        
        %%%% range on x vector
        xsub1 = find(x >= subr1(1) & x <= subr1(2)); 
        xsub2 = find(x >= 0 & x <= subr2(2));     
        
        %%%% indices on x vector
        Isub1_l = xsub1(1)-2 ; Isub1_r = size(x,2); 
        Isub2_l = 1 ; Isub2_r = xsub2(end)+2; 
        
        %%% Taking -2 is safe because the lower end is not close to edge.        
        %%%%  the two subregions around centrosome 1. 
        a_sub1 = th(th >= ang1(1));
        a_sub2 = th(th <= ang1(2)); 
        
        %%%% overlap region along s 
        s_sub1 = round(s(th >= ang1(1)),rs);  
        s_sub2 = s(th <= ang1(2));
        
        %%%% correct for the possibility that s_sub1(end) or last few
        %%%% points may be > L due to discretisation. e.g.  L =20.2 while
        %%%% s_sub1(end) = 20.23. This will give interpolation errors.
        %%%% This possibility does not arise in s_sub2 because the lowest
        %%%% point is 0 in both x and s grids. 
        s_sub1 = s_sub1(s_sub1<= sys_L);  
       
        %%%% interpolate concentration on s
        c_sub1 = interp1(x(Isub1_l:Isub1_r),conc(Isub1_l:Isub1_r),s_sub1); 
        c_sub2 = interp1(x(Isub2_l:Isub2_r),conc(Isub2_l:Isub2_r),s_sub2); 
   
        %%%% checks if the interpolation is done right
        chkp1 = find(isnan(c_sub1));
        chkp2 = find(isnan(c_sub2));
        
        if (size(chkp1,2)>=1 || size(chkp2,2)>=1 )
            error('Interpolation error in 2 crossings, centrosome 1, case 2')
        end
    
        den_sub1 = sqrt((r(a_sub1).*cos(a_sub1)-lsp*cos(phi1)).^2+...
                   (r(a_sub1).*sin(a_sub1)-lsp*sin(phi1)).^2);
               
        torq_sub1 = lsp*r(a_sub1).*c_sub1.*sin(a_sub1-phi1)./den_sub1;
    
        den_sub2 = sqrt((r(a_sub2).*cos(a_sub2)-lsp*cos(phi1)).^2+...
                   (r(a_sub2).*sin(a_sub2)-lsp*sin(phi1)).^2);
               
        torq_sub2 = lsp*r(a_sub2).*c_sub2.*sin(a_sub2-phi1)./den_sub2;
   
        torq_dens1(th >= ang1(1)) = torq_sub1; 
        torq_dens1(th <= ang1(2)) = torq_sub2; 
        
        torq1 = trapz(torq_sub1*dth) + trapz(torq_sub2*dth); 
   
    end
       
    %%%%% Calculate the torque on centrosome 2 similarly. There is a
    %%%%% possiblity of positive region crossing the boundary.     
    ang2 = round([th(Ivec2(1)+1) th(Ivec2(end))],rn);    

    if fun2(th(Ivec2(1)+2)) > 0  
        Exch2 =0;
    else
        Exch2 =1;
    end
    
    if Exch2 ==1
    [ang2] = [ang2(2) ang2(1)]; 
    end
    
    
    %%%%%%%%%%%% CASE 1 :  NO CROSSING OF BOUNDARIES %%%%%%%%%%%%%%%%%%%%%%
    if (Exch2 == 0)  
        s2 = [s(Ivec2(1)+1) s(Ivec2(end))];
    
        %%% input range for interpolation 
        xran2 = find(x >= s2(1) & x <= s2(2)); 
        Il2 = xran2(1)-1 ; Ir2 = xran2(end)+1; 
        
        %%%% specification in case xran2 end points are exactly on one of
        %%%% the edges of linear system.
        if Il2 ==0;  Il2 = 1; end
        if Ir2 == size(x,2)+1; Ir2 = size(x,2); end
    
        %%% output range for interpolation 
        arange2 = th(th > ang2(1) & th < ang2(2));
        srange2 = s(th > ang2(1) & th < ang2(2));
    
        %%%% interpolate concentration on angular/peripheral grid
        c2interp = interp1(x(Il2:Ir2),conc(Il2:Ir2),srange2); 
    
        chk2 = find(isnan(c2interp));
        if (size(chk2,2)>=1 )
            error('Interpolation error in 2 crossings, centrosome 2, case 1')
        end
     
        den2 = sqrt((r(arange2).*cos(arange2)-lsp*cos(phi2)).^2+...
                   (r(arange2).*sin(arange2)-lsp*sin(phi2)).^2);
      
        torq_dens2(th > ang2(1) & th < ang2(2)) = lsp*r(arange2).*c2interp.*sin(arange2-phi2)./den2;  
        torq2 = trapz(dth*lsp*r(arange2).*c2interp.*sin(arange2-phi2)./den2); 
    
     elseif(Exch2 == 1)
        
        subr1 = [s(Ivec2(end)) x(end)];    %% region in 4th quad 
        subr2 = [0 s(Ivec2(1)+1)];         %% region 1st quad        
        
        %%%% range on x vector
        xsub1 = find(x >= subr1(1) & x <= subr1(2)); 
        xsub2 = find(x >= 0 & x <= subr2(2));     
        
        %%%% indices on x vector
        Isub1_l = xsub1(1)-2 ; Isub1_r = size(x,2); 
        Isub2_l = 1 ; Isub2_r = xsub2(end)+2; 
        
        %%% Taking -2 is safe because the lower end is not close to edge.        
        %%%%  the two subregions around centrosome 1. 
        a_sub1 = th(th >= ang2(1));
        a_sub2 = th(th <= ang2(2)); 
        
        %%%% overlap region along s 
        s_sub1 = round(s(th >= ang2(1)),rs);  
        s_sub2 = s(th <= ang2(2));
        
        %%%% correct for the possibility that s_sub1(end) or last few
        %%%% points may be > L due to discretisation. e.g.  L =20.2 while
        %%%% s_sub1(end) = 20.23. This will give interpolation errors.
        %%%% This possibility does not arise in s_sub2 because the lowest
        %%%% point is 0 in both x and s grids. 
        s_sub1 = s_sub1(s_sub1<= sys_L);  
       
        %%%% interpolate concentration on s
        c_sub1 = interp1(x(Isub1_l:Isub1_r),conc(Isub1_l:Isub1_r),s_sub1); 
        c_sub2 = interp1(x(Isub2_l:Isub2_r),conc(Isub2_l:Isub2_r),s_sub2); 
   
        %%%% checks if the interpolation is done right
        chkp1 = find(isnan(c_sub1));
        chkp2 = find(isnan(c_sub2));
        
        if (size(chkp1,2)>=1 || size(chkp2,2)>=1 )
            error('Interpolation error in 2 crossings, centrosome 2, case 2')
        end
    
        den_sub1 = sqrt((r(a_sub1).*cos(a_sub1)-lsp*cos(phi2)).^2+...
                   (r(a_sub1).*sin(a_sub1)-lsp*sin(phi2)).^2);
               
        torq_sub1 = lsp*r(a_sub1).*c_sub1.*sin(a_sub1-phi2)./den_sub1;
    
        den_sub2 = sqrt((r(a_sub2).*cos(a_sub2)-lsp*cos(phi2)).^2+...
                   (r(a_sub2).*sin(a_sub2)-lsp*sin(phi2)).^2);
               
        torq_sub2 = lsp*r(a_sub2).*c_sub2.*sin(a_sub2-phi2)./den_sub2;
   
        torq_dens1(th >= ang2(1)) = torq_sub1; 
        torq_dens1(th <= ang2(2)) = torq_sub2; 
        
        torq2 = trapz(torq_sub1*dth) + trapz(torq_sub2*dth); 
     end
    
    
    %%%%%%  Calculate resulting torque on the spindle. 
    Omega = torq1 + torq2;
    
    %%%% plot if flag is passed.
    if flag == true
    reg1 = (th >= ang1(1) & th < ang1(2)); 
    reg2 = (th >= ang2(1) & th < ang2(2));  
    
    if (Exch1 == 1); reg1 = ( th <= ang1(2) || th >= ang1(1)); end 
    
    if (Exch2 == 2); reg2 = ( th <= ang2(2)|| th >= ang2(1)); end 
       
    figure(1) 
    hold on; plot(xunit(reg1),yunit(reg1),'bo','Markersize',3);
    hold on; plot(xunit(reg2),yunit(reg2),'bo','Markersize',3);
    axis equal tight          
    end
                   
%%%%%%%% 4 CROSSING BEGINS HERE : TWO CASES POSSIBLE FOR 4 CROSSINGS TOO %%  
%%%%%%%% The 4  OVERLAP REGIONS ARE IN 4 DIFFERENT QUADS %%%%%%%%%%%%%%%%%%
%%%%%%%% Like 2 crossings case, the crossing of linear boundaries possible
%%%%%%%% here too. In addition this can be done by region 1 or region 2
%%%%%%%% around the centrosome 1 : corresponds to crossings from left or right. 
%%%% Sometimes around centrosome 1, wrong pairs are picked up as the two regions. 
%%%% This is a mistake and can be corrected by changing the order of ang1 
%%%% entries. This is not a problem for centrosome 2
%%%% because of the handedness in the definition of theta grid. The
%%%% first two angles of ang2 always corresond to correct region.

elseif (size(Ivec1,2)>2)   
    
    ang1 = round([th(Ivec1(1)+1) th(Ivec1(2))  th(Ivec1(3)+1)  th(Ivec1(4))],rn);
    s1 = [s(Ivec1(1)+1) s(Ivec1(2))  s(Ivec1(3)+1)  s(Ivec1(4))] ;
    
    %%% condition that region 1 crosses left linear boundary and is in
    %%% first and 4 angular quad :  dealt with in region 1, case 2
    
    cond1 = (ang1(1) >=0 && ang1(1) < pi) && ...
           (ang1(2) >=pi && ang1(2)<= 2*pi) && ... 
           (ang1(3) >=pi && ang1(3)<= 2*pi) && ...
           (ang1(4) >=pi && ang1(4)<= 2*pi); 
   
    %%% condition that region 2 crosses right linear boundary and is in
    %%% first and 4 angular quad :  dealt with in region 2, case 2
    
    cond2 = (ang1(1) >=0 && ang1(1) < pi) && ...
           (ang1(2) >=0 && ang1(2)<= pi) && ... 
           (ang1(3) >=0 && ang1(3)<= pi) && ...
           (ang1(4) >=pi && ang1(4)<=2*pi); 

%     %%%% correct for the order of angles
     if cond1 
         ang1 = [ang1(4) ang1(1) ang1(2) ang1(3)];
         s1 = [s1(4) s1(1) s1(2) s1(3)]; 
     end
     
     if cond2
          ang1 = [ang1(2) ang1(3) ang1(4) ang1(1)];
          s1 = [s1(2) s1(3) s1(4) s1(1)]; 
     end
    
    %%%% torque on centrosome 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (sin(ang1(1))*sin(ang1(2)) > 0  )
    
        %%%% range on x vector
        xran11 = find(x >= s1(1) & x <= s1(2)); 
   
        %%%% indices on x vector:  input for interpolation
        Il11 = xran11(1)-1 ; Ir11 = xran11(end)+1; 

        if xran11(1) ==1
            Il11 =1;
        end
        %%%%  angular overlap regions around centrosome 1 and 2
        arange11 = th(th >= ang1(1) & th <= ang1(2));

        %%%% overlap region along s: output for interpolation
        srange11 = s(th >= ang1(1) & th <= ang1(2));
 
        %%%% interpolate concentration on s 
        c11interp = interp1(x(Il11:Ir11),conc(Il11:Ir11),srange11); 
    
        %%%% checks if the interpolation is done right
        chk11 = find(isnan(c11interp));
  
        if (size(chk11,2)>=1 )
            error('Interpolation error: 4 crossings,region1,case 1')
        end

       %%% calculate torques from each region.
       den11 = sqrt((r(arange11).*cos(arange11)-lsp*cos(phi1)).^2+...
                   (r(arange11).*sin(arange11)-lsp*sin(phi1)).^2);
               
       torq_dens1(th >= ang1(1) & th <= ang1(2))= lsp*r(arange11).*c11interp.*sin(arange11-phi1)./den11;
       
       torq11 = trapz(dth*lsp*r(arange11).*c11interp.*sin(arange11-phi1)./den11);
   
    
    elseif (sin(ang1(1))*sin(ang1(2)) <= 0  )
 
        sub11_lo = [s(Ivec1(4)) x(end)] ;     %%% branch in quad 4
        sub11_hi = [0 s(Ivec1(1)+1)] ;        %%% branch in quad 1
       
        %%%% range on x vector
        xsub11_lo = find(x >= sub11_lo(1) & x <= sub11_lo(2)); 
        xsub11_hi = find(x >= 0 & x <= sub11_hi(2)); 
    
        %%%% indices on x vector
        Isub11_lo_l = xsub11_lo(1)-1 ; Isub11_lo_r = size(x,2); 
        Isub11_hi_l = 1 ;              Isub11_hi_r = xsub11_hi(end)+1; 
  
        %%%%  the two subregions around centrosome 1. 
        a_sub_11_lo = th(th >= ang1(1));
        a_sub_11_hi = th(th <= ang1(2));
 
        %%%% overlap region along s 
        s_sub_11_lo = round(s(th >= ang1(1)),rs);
        s_sub_11_hi = s(th <= ang1(2));
        
        if s_sub_11_lo(end)>sys_L
            s_sub_11_lo = s_sub_11_lo(1:end-1);  
        end
        
        
        %%%% interpolate concentration on s 
        c_sub_11_lo = interp1(x(Isub11_lo_l:Isub11_lo_r),conc(Isub11_lo_l:Isub11_lo_r),s_sub_11_lo); 
        c_sub_11_hi = interp1(x(Isub11_hi_l:Isub11_hi_r),conc(Isub11_hi_l:Isub11_hi_r),s_sub_11_hi); 
   
        %%%% checks if the interpolation is done right
        chkp1 = find(isnan(c_sub_11_lo));
        chkp2 = find(isnan(c_sub_11_hi));
        
        if (size(chkp1,2)>=1 || size(chkp2,2)>=1 )
            error('Interpolation error: 4 crossings, region 1, case 2')
        end
        
        den_sub_11_lo = sqrt((r(a_sub_11_lo).*cos(a_sub_11_lo)-lsp*cos(phi1)).^2+...
                   (r(a_sub_11_lo).*sin(a_sub_11_lo)-lsp*sin(phi1)).^2);
               
        torq_sub_11_lo = trapz(dth*lsp*r(a_sub_11_lo).*c_sub_11_lo.*sin(a_sub_11_lo - phi1)./den_sub_11_lo);
        
        den_sub_11_hi = sqrt((r(a_sub_11_hi).*cos(a_sub_11_hi)-lsp*cos(phi1)).^2+...
                   (r(a_sub_11_hi).*sin(a_sub_11_hi)-lsp*sin(phi1)).^2);
               
        torq_sub_11_hi = trapz(dth*lsp*r(a_sub_11_hi).*c_sub_11_hi.*sin(a_sub_11_hi - phi1)./den_sub_11_hi);
       
        torq_dens1(th >= ang1(1)) = lsp*r(a_sub_11_lo).*c_sub_11_lo.*sin(a_sub_11_lo - phi1)./den_sub_11_lo; 
        torq_dens1(th <= ang1(2)) = lsp*r(a_sub_11_hi).*c_sub_11_hi.*sin(a_sub_11_hi - phi1)./den_sub_11_hi; 
        
        torq11 = torq_sub_11_hi + torq_sub_11_lo; 
        
    end
     
    %%%%  torque on centrosome 1 from region 2. 
    if (sin(ang1(3))*sin(ang1(4)) > 0 )   
    
        %%%% range on x vector
        xran12 = find(x >= s1(3) & x <= s1(4)); 
        
        %%%% indices on x vector:  input for interpolation
        Il12 = xran12(1)-1 ; Ir12 = xran12(end)+1; 
        
        %%%%  angular overlap regions around centrosome 1 and 2
        arange12 = th(th >= ang1(3) & th <= ang1(4));

        %%%% overlap region along s: output for interpolation
        srange12 = s(th >= ang1(3) & th <= ang1(4));
 
        %%%% interpolate concentration on s 
        c12interp = interp1(x(Il12:Ir12),conc(Il12:Ir12),srange12);
    
        %%%% checks if the interpolation is done right
        chk12 = find(isnan(c12interp));
  
        if (size(chk12,2)>=1 )
            error('Interpolation error: 4 crossings,region2,case 1')
        end

        den12 = sqrt((r(arange12).*cos(arange12)-lsp*cos(phi1)).^2+...
                   (r(arange12).*sin(arange12)-lsp*sin(phi1)).^2);
               
        torq_dens1(th >= ang1(3) & th <= ang1(4)) = lsp*r(arange12).*c12interp.*sin(arange12-phi1)./den12; 
               
        torq12 = trapz(dth*lsp*r(arange12).*c12interp.*sin(arange12-phi1)./den12);
    
    
    elseif (sin(ang1(3))*sin(ang1(4)) <= 0 )
       
        sub_12_lo = [s(Ivec1(3)+1) x(end)] ;   %%% branch in quad 4
        sub_12_hi = [0 s(Ivec1(4))] ;        %%% branch in quad 1
       
        %%%% range on x vector
        xsub_12_lo = find(x >= sub_12_lo(1) & x <= sub_12_lo(2)); 
        xsub_12_hi = find(x >= 0 & x <= sub_12_hi(2)); 
    
        %%%% indices on x vector
        Isub12_lo_l = xsub_12_lo(1)-1 ; Isub12_lo_r = size(x,2); 
        Isub12_hi_l = 1 ;               Isub12_hi_r = xsub_12_hi(end)+1; 
  
        %%%%  the two subregions around centrosome 1. 
        a_sub_12_lo = th(th >= ang1(3));
        a_sub_12_hi = th(th <= ang1(4));
 
        %%%% overlap region along s 
        s_sub_12_lo = round(s(th >= ang1(3)),rs);
        s_sub_12_hi = s(th <= ang1(4));
        
        if s_sub_12_lo(end)>sys_L
            s_sub_12_lo = s_sub_12_lo(1:end-1);  
        end
   
        %%%% interpolate concentration on s 
        c_sub_12_lo = interp1(x(Isub12_lo_l:Isub12_lo_r),conc(Isub12_lo_l:Isub12_lo_r),s_sub_12_lo); 
        c_sub_12_hi = interp1(x(Isub12_hi_l:Isub12_hi_r),conc(Isub12_hi_l:Isub12_hi_r),s_sub_12_hi); 
   
        %%%% checks if the interpolation is done right
        chkp1 = find(isnan(c_sub_12_lo));
        chkp2 = find(isnan(c_sub_12_hi));
        
        if (size(chkp1,2)>=1 || size(chkp2,2)>=1 )
            error('Interpolation error: 4 crossings, region 2, case 2')
        end
    
        den_sub_12_lo = sqrt((r(a_sub_12_lo).*cos(a_sub_12_lo)-lsp*cos(phi1)).^2+...
                   (r(a_sub_12_lo).*sin(a_sub_12_lo)-lsp*sin(phi1)).^2);
               
        torq_sub_12_lo = trapz(dth*lsp*r(a_sub_12_lo).*c_sub_12_lo.*sin(a_sub_12_lo-phi1)./den_sub_12_lo);
        
        den_sub_12_hi = sqrt((r(a_sub_12_hi).*cos(a_sub_12_hi)-lsp*cos(phi1)).^2+...
                   (r(a_sub_12_hi).*sin(a_sub_12_hi)-lsp*sin(phi1)).^2);
               
        torq_sub_12_hi = trapz(dth*lsp*r(a_sub_12_hi).*c_sub_12_hi.*sin(a_sub_12_hi-phi1)./den_sub_12_hi);
        
        torq_dens1(th >= ang1(3)) = lsp*r(a_sub_12_lo).*c_sub_12_lo.*sin(a_sub_12_lo-phi1)./den_sub_12_lo; 
        torq_dens1(th <= ang1(4)) = lsp*r(a_sub_12_hi).*c_sub_12_hi.*sin(a_sub_12_hi-phi1)./den_sub_12_hi; 
        
        torq12 = torq_sub_12_lo + torq_sub_12_hi; 
        
    end
        
        %%%% Torque on centrosome 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        
        ang2 = round([th(Ivec2(1)+1) th(Ivec2(2))  th(Ivec2(3)+1)  th(Ivec2(4))],rn) ;
        s2 = [s(Ivec2(1)+1) s(Ivec2(2))  s(Ivec2(3)+1)  s(Ivec2(4))] ;   
        xran21 = find(x >= s2(1) & x <= s2(2)); 
        xran22 = find(x >= s2(3) & x <= s2(4)); 
        
        %%%% indices on x vector:  input for interpolation
        Il21 = xran21(1)-1 ; Ir21 = xran21(end)+1; 
        Il22 = xran22(1)-1 ; Ir22 = xran22(end)+1;
    
        %%%%  angular overlap regions around centrosome 1 and 2
        arange21 = th(th >= ang2(1) & th <= ang2(2));
        arange22 = th(th >= ang2(3) & th <= ang2(4));
    
        %%%% overlap region along s: output for interpolation
        srange21 = s(th >= ang2(1) & th <= ang2(2));
        srange22 = s(th >= ang2(3) & th <= ang2(4));

        %%%% interpolate concentration on s 
        c21interp = interp1(x(Il21:Ir21),conc(Il21:Ir21),srange21); 
        c22interp = interp1(x(Il22:Ir22),conc(Il22:Ir22),srange22); 
    
        %%%% checks if the interpolation is done right
    
        chk21 = find(isnan(c21interp));
        chk22 = find(isnan(c22interp));
    
        if (size(chk21,2)>=1 || size(chk22,2)>=1 )
        error('Interpolation error: check 4 crossings, centrosome 2')
        end
        
        den21 = sqrt((r(arange21).*cos(arange21) - lsp*cos(phi2)).^2+...
                                   (r(arange21).*sin(arange21) - lsp*sin(phi2)).^2);
               
        torq_dens2(th >= ang2(1) & th <= ang2(2)) = ...
                            lsp*r(arange21).*c21interp.*sin(arange21 - phi2)./den21; 
               
        torq21 = trapz(dth*lsp*r(arange21).*c21interp.*sin(arange21 - phi2)./den21);
    
        den22 = sqrt((r(arange22).*cos(arange22) - lsp*cos(phi2)).^2+...
                                   (r(arange22).*sin(arange22) - lsp*sin(phi2)).^2);
        
        torq_dens2(th >= ang2(3) & th <= ang2(4)) = ...
                            lsp*r(arange22).*c22interp.*sin(arange22 - phi2)./den22; 
                    
        torq22 = trapz(dth*lsp*r(arange22).*c22interp.*sin(arange22 - phi2)./den22);
   
 
    Omega = torq11+torq12+torq21+torq22;
    
    
    if flag == true 
    reg1 = (th >= ang1(1) & th < ang1(2)); 
    if (sin(ang1(1))*sin(ang1(2))) < 0
        reg1 = (th >= ang1(1) | th < ang1(2));
    end
    
    reg2 = (th >= ang1(3) & th < ang1(4));  
    if (sin(ang1(3))*sin(ang1(4))) < 0
        reg2 = (th >= ang1(3) | th < ang1(4));
    end
    
    reg3 = (th >= ang2(1) & th < ang2(2));  
    reg4 = (th >= ang2(3) & th < ang2(4));  
    
    hold on; plot(xunit(reg1),yunit(reg1),'bo',xunit(reg2),yunit(reg2),'bo','Markersize',3)
    hold on; plot(xunit(reg3),yunit(reg3),'bo',...
          xunit(reg4),yunit(reg4),'bo','Markersize',3);
    drawnow
 
    axis equal tight
    end
    
    
end
% torq_dens = torq_dens1 + torq_dens2; 







































