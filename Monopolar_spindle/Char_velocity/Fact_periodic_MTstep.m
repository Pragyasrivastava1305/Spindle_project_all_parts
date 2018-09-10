function [rhs, Ipos ,Ineg] = Fact_periodic_MTstep(xcen,Lm,x,dx,dx_nuc,conc)
%%%% given the profile of LGN this function calculates the force on
%%%% centrosome with periodic boundary conditions. periodicity needs to be
%%%% applied to MT branches. To calculate force, conc profile is
%%%% interpolated on a grid in MT region about the Centrosome. 
%%%% The range of interpolation on initial grid should be a bit bigger than
%%%% the range on final grid for interpolation to work. 
%%%% The code allows the grid-spacing to be different on MT branches from
%%%% that of the spatial grid on which concentration is defined. These two
%%%% are denoted by dx_nuc and dx. 

%%%% The algorithm is to first check if both MT branches are in same
%%%% periodic image. If not, then the branches are wrapped back into same
%%%% periodic image from left or right, depending on which branch is
%%%% outside the periodic image. Then concentration is interpolated on this
%%%% interpolated and wrapped grid. In final stape, the integration is
%%%% performed on right and left branches to calculated positive and
%%%% negative forces. 

Nend=size(x,2);     ninter=1+2*Lm/dx_nuc;

%%%%% centrosome position and interpolation range on final grid
xcen_rnd=round(xcen,2); 
xlo=xcen_rnd-Lm; xhi=xcen_rnd+Lm; 

%%%% nearest index to centrosome position on initial grid 
xcen_ind=1+floor(xcen_rnd/dx); 



%%%% Case 1: both branches in same periodic image

if (xlo>=0 && xhi <= x(end))    

%%%% Range of final grid.  It is chosen slightly less than the MT range 
%%%% first point is xlo-dx_nuc and xhi+dx_nuc so that this length is always
%%%% in interpolation range.  There will be inaccuracy in force calculation
%%%% over the length scale of 2*dx_nuc.
    
xfin=xlo+dx_nuc:dx_nuc:xhi-dx_nuc;    


Ilo=floor(xlo/dx)+1;           %%% Floor and ceil functions will capture 
Ihi=1+ceil(xhi/dx);            %%% the requirement of initial range being 
                               %%% bigger than final range. 

xinit=x(Ilo:Ihi);                             
cinit=conc(Ilo:Ihi);

cfin_p=interp1(xinit,cinit,xfin); 
% ninter=size(xfin,2); 


Ipos=xcen_ind:Ihi;             %%%% extent of positive arm on diff grid
Ineg=Ilo:xcen_ind;             %%%% extent of negative arm on diff grid
                               %%%% Only used for plotting and not for
                               %%%% force calculation

fpos=dx_nuc*trapz(cfin_p(1+(ninter-1)/2:end));  %%% force calculation on  
fneg=dx_nuc*trapz(cfin_p(1:1+(ninter-1)/2));    %%% interpolated grid.


%%%% right branch in next periodic image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% force calculation on left branch remains same. Force calculation on
%%%% right breaks in two parts. 

elseif (xlo>=0 && xhi>x(end))   
%%%% separating the parts in two different periodic images for proper
%%%% interpolation. 
   
xfin_p1=xlo:dx_nuc:x(end)-dx_nuc;               %%%% in first image
xfin_p2=x(1):dx_nuc:xhi-x(end);                 %%%% in second image

% ninter=size(xfin_p1,2)+size(xfin_p2,2); 
% ninter=3001;

Ilo1=1+floor(xfin_p1(1)/dx);            Ihi1=1+ceil(xfin_p1(end)/dx);
xinit1=x(Ilo1-1:Ihi1);                   
cinit1=conc(Ilo1-1:Ihi1);

%%%%(a) an extra -1 to ensure that the original grid is bigger than the
%%%% interpolating grid. 


Ilo2=1+floor(xfin_p2(1)/dx);            Ihi2=1+ceil(xfin_p2(end)/dx);

if (size(xfin_p2,2)==1)   %%%% this condition takes care of the case when 
Ihi2=2;                   %%%% xcen+Lm-L=0.01  and there is only 1 point in 
end                       %%%% next periodic image. 

%%% Reason for an extra +1 same as (a) on the other branch. 

xinit2=x(Ilo2:Ihi2+1);
cinit2=conc(Ilo2:Ihi2+1);


cfin_p1=interp1(xinit1,cinit1,xfin_p1); 
cfin_p2=interp1(xinit2,cinit2,xfin_p2); 

%%%%% For plotting 
    Ipos1=xcen_ind:Nend; 
    Ipos2=1:Ihi2;
    Ineg=Ilo1:xcen_ind;
    
    Ipos=zeros(1,size(Ipos1,2)+size(Ipos2,2)); 
    Ipos(1,1:size(Ipos1,2))=Ipos1; Ipos(1,size(Ipos1,2)+1:end)=Ipos2;


%%%%% For force calculation
fneg=dx_nuc*trapz(cfin_p1(1:(ninter-1)/2));   %%% fneg remains same

%%%% part in first periodic image
Iind=1+(ninter-1)/2:size(xfin_p1,2);

fpos1=dx_nuc*trapz(cfin_p1(Iind)); 

%%%% part in second periodic image
fpos2=dx_nuc*trapz(cfin_p2);  

fpos=fpos1+fpos2; 

%%%% left branch in previous periodic image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif (xlo<0 && xhi<=x(end))   

xfin_p1=x(1):dx_nuc:xhi;                               %%%% in first image
xfin_p2=x(end)+xlo:dx_nuc:x(end)-dx_nuc;               %%%% in previous image
% ninter=size(xfin_p1,2)+size(xfin_p2,2); 
% ninter=3001;

Ilo1=1+floor(xfin_p1(1)/dx);            Ihi1=1+ceil(xfin_p1(end)/dx);
xinit1=x(Ilo1:Ihi1+1);                   
cinit1=conc(Ilo1:Ihi1+1);

Ilo2=1+floor(xfin_p2(1)/dx);            Ihi2=1+ceil(xfin_p2(end)/dx);
xinit2=x(Ilo2-1:Ihi2);
cinit2=conc(Ilo2-1:Ihi2);

cfin_p1=interp1(xinit1,cinit1,xfin_p1); 
cfin_p2=interp1(xinit2,cinit2,xfin_p2); 


%%%%% for plotting 
    Ipos=xcen_ind:Ihi1; 
    Ineg1=1:xcen_ind;
    Ineg2=Ilo2:Ihi2;
    
    Ineg=zeros(1,size(Ineg1,2)+size(Ineg2,2)); 
    Ineg(1,1:size(Ineg1,2))=Ineg1; Ineg(1,size(Ineg1,2)+1:end)=Ineg2;


%%%%% For force calculation
fpos=dx_nuc*trapz(cfin_p1(end-(ninter-1)/2:end)); %%% fpos remains same

%%%% part in first periodic image
Jind=1:size(xfin_p1,2)-(ninter-1)/2;

fneg1=dx_nuc*trapz(cfin_p1(Jind)); 

%%%% part in second periodic image
fneg2=dx_nuc*trapz(cfin_p2);  

fneg=fneg1+fneg2; 

end

rhs=fpos-fneg; 


end

