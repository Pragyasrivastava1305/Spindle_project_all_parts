function [rhs,IR,IL,xwrap,Pdfwrap,cs] = Fact_periodic_Mtdist(xcen,xMT,PdfMT,dx_MT,x,dx,r,conc)
%%%% This is improved version of the force calculation using MT end
%%%% distribution than the one used for experimental data. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% xcen = current position of centrosome
%%% xMT = xpoints of the Pdf of MTs. 
%%% PdfMT = Pdf of Mt ends 
%%% x = spatial grid on which conc is defined
%%% dx = gridsize of LGN grid. 
%%% conc profile of LGN as given by experimental or simulation data. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IR = region that applies force in positive direction
%%% IL =  region that applies force in negative direction
%%% rhs = Active pulling force.
%%% xnew = new x range (wrapped) for MT distribution
%%% Pdfnew = new Pdf of MTs (wrapped) for MT distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% From xcen and dx_nuc construct a symmetric region around centrosome
%%%% where MTs are present. Interpolate PdfMT on this region. Interpolate
%%%% concentration on this region. Determine positive and negative
%%%% branches. Integrate to calculate forces. When MT branches are not
%%%% contained in same periodic image, a wrapping-back to the same image is
%%%% performed. Depending on whether the left or right arm crosses the
%%%% system boundaries, all 5 cases are considered : 1. both arms in same
%%%% periodic image, 2. right/left arm exactly on the boundary, 3.
%%%% right/left arm in next/previous periodic image. Some lines of this
%%%% code give no error in force calculation when the force is assumed to
%%%% come from the MT ends. If the force is assumed is to come from all
%%%% over the length of MTs then the same lines will give error in forces
%%%% over a distance of dx. This may be significant and may be strong
%%%% enough to induce an asymmetry.  

% xcen = xC_per_new; 
% xMT = MT_x;  PdfMT = MT_pdf; dx_MT = dx_nuc; 
% conc = cL_new; 

if dx_MT > dx
    warning('dx_nuc >= dx, may lead to some variables being out of bounds')
end


%%%%%% CREATE SYMMETRIC PROFILE OF MTs AROUND THE CENTROSOME %%%%%%%%%%%%%%%
PdfMT = PdfMT/sum(PdfMT);  

%%%%% round off centrosome position and define the grid around centrosome 
%%%%% with resolution dx_MT. 
xcen_rnd=round(xcen,r);  
xcen_ind=1+floor(xcen_rnd/dx);       %%%%%%%%% index of xcen on LGN grid

xMT = round(xMT,r);

xlo= xcen_rnd - xMT(end);        
xhi= xcen_rnd + xMT(end); 
xhires = xMT(1):dx_MT:xMT(end);      %%%% increase resolution of PDF x-axis
Pdf_interp = interp1(xMT,PdfMT,xhires);          %%%% PDF defined on xhires

%%% Define symmetric region and MT distribution around centrosome.
 xnew = xlo:dx_MT:xhi;   

%%%%% index of centrosome position on MT grid
if mod(size(xnew,2),2)~=0                
    Ncen=(1+size(xnew,2))/2;     
else 
    Ncen=size(xnew,2)/2; 
end

Pdfnew = zeros(1,size(xnew,2)); 
%%%%% Create symmetric profile of MTs around centrosome %%%%%%%%%%%%%%%%%%%
Pdfnew(end-size(Pdf_interp,2)+1:end) = Pdf_interp;
Pdfnew(1:size(Pdf_interp,2)) = fliplr(Pdf_interp); 

%%%%% This Pdf has to be wrapped back in the same periodic image. Following
%%%%% lines of code perform this, interpolate the concentration field on
%%%%% wrapped up arms and then integrate it using PDF of MTs as weight. 

%%%%%  PERIODIC WRAPPING AND FORCE CALCULATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% We need to consider three separate cases.  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% CASE 1: both branches in same periodic image %%%%%%%%%%%%%%%%%%%%%%%%%

if (xlo>=dx && xhi <= x(end)-dx)    
    cs = 'case 1'; 
    
    Ilo_grid1 = floor(xlo/dx)+1;             %%% 1 added as floor(xlo/dx)=0
    Ihi_grid1 = 1+ceil(xhi/dx);

    IL = Ilo_grid1:xcen_ind ; 
    IR = xcen_ind:Ihi_grid1 ; 

    cL= interp1(x(IL(1)-1:IL(end)+1),conc(IL(1)-1:IL(end)+1),xnew(1:Ncen)); 
    cR= interp1(x(IR(1)-1:IR(end)+1),conc(IR(1)-1:IR(end)+1),xnew(Ncen+1:end)); 
    
    fneg=dx_MT*trapz(cL(1:Ncen-1).*Pdfnew(1:Ncen-1));    
    fpos=dx_MT*trapz(cR.*Pdfnew(Ncen+1:end));  

    xwrap = xnew; 
    Pdfwrap = Pdfnew;                       %%% interpolated grid.

%%%% CASE 2 : end of right branch exactly on L. only 1 point in next
%%%% periodic image. There is no need to wrap but care must be taken to
%%%% interpolate on the last point. The value of concentration on last
%%%% point of xnew then is taken to be exactly the concentration at x(end).
elseif (xhi > x(end)-dx && xhi <= x(end))

    cs = 'case 2';
    Ilo_grid1 = floor(xlo/dx)+1; 
    Ihi_grid1 = 1+ceil(xhi/dx);

    IL = Ilo_grid1:xcen_ind ; 
    IR = xcen_ind:Ihi_grid1; 

    cL= interp1(x(IL(1)-1:IL(end)+1),conc(IL(1)-1:IL(end)+1),xnew(1:Ncen));
    cR= zeros(1,size(xnew(Ncen+1:end),2));
    cR(1:end-1)= interp1(x(IR(1)-1:IR(end)),conc(IR(1)-1:IR(end)),xnew(Ncen+1:end-1));
    cR(end) = conc(IR(end)); 

    fneg=dx_MT*trapz(cL(1:Ncen-1).*Pdfnew(1:Ncen-1));    
    fpos=dx_MT*trapz(cR.*Pdfnew(Ncen+1:end));  

    xwrap = xnew; 
    Pdfwrap = Pdfnew;   

%%%% CASE 3: end of left branch exactly at the boundary. Same as case 2 but
%%%% on left boundary

elseif (xlo>=x(1) && xlo < x(1)+dx)

    cs = 'case 3';

    Ilo_grid1 = floor(xlo/dx)+1;                %%% 1 added as floor(xlo/dx)=0. 
    Ihi_grid1 = 1+ceil(xhi/dx);

    IL = Ilo_grid1:xcen_ind ; 
    IR = xcen_ind:Ihi_grid1; 

    cL= zeros(1,size(xnew(1:Ncen),2));
    cL(2:end)= interp1(x(IL(1):IL(end)+1),conc(IL(1):IL(end)+1),xnew(2:Ncen));
    cL(1)= conc(IL(1)); 
    cR= interp1(x(IR(1)-1:IR(end)),conc(IR(1)-1:IR(end)),xnew(Ncen+1:end));

    fneg=dx_MT*trapz(cL(1:Ncen-1).*Pdfnew(1:Ncen-1));    
    fpos=dx_MT*trapz(cR.*Pdfnew(Ncen+1:end));  

    xwrap = xnew; 
    Pdfwrap = Pdfnew;   

%%%%% CASE 4 : right branch in next periodic image. left branch remains unchanged. 
elseif (xlo>x(1) && xhi>x(end))   
    cs = 'case 4'; 
    %%%%%  right branch needs to be wrapped. 
    %%%%%  need to redefine xnew and Pdfnew. Part of both these will be in
    %%%%%  right periodic image. To do this, define xwrap and pdfwrap. 
    lenx = xnew(end) - x(end); 
    Ip = find(xnew >= x(end));   %%% find indices of the part of xnew that is in next image

    xwrap = zeros(1,size(xnew,2));
    Pdfwrap = zeros(1,size(xnew,2)); 
    xwrap(1:size(Ip,2)) = xnew(Ip)-x(end);
    xwrap(size(Ip,2)+1:end) = xnew(1:Ip(1)-1);

    %%%%% Ip(1), because xnew(Ip(1)-1)=x(end) but this is already counted
    %%%%% as first point of xwrap in previous line. x(1) =x(end) 
    Pdfwrap(1:size(Ip,2)) = Pdfnew(Ip);
    Pdfwrap(size(Ip,2)+1:end) = Pdfnew(1:Ip(1)-1);

    %%%%%%% force calculation from left branch remains same as case 1  
    %%% on LGN grid
    Ilo_grid1 = floor(xlo/dx)+1; 
    IL = Ilo_grid1:xcen_ind ;   
 
    %%%% left branch starts at xlo which is the (size(Ip,2)+1)th point 
    %%%% of xwrap. and it goes untill Ncen length on xwrap grid.  When
    %%%% centrosome is exactly at the edge i.e. xcen = x(end), then IL(end)+1
    %%%% > x(end). This gives trouble so we treat this case separately. 

    Left_Ind = size(Ip,2)+1:size(Ip,2)+Ncen-1; 

    if (x(IL(end))~=x(end)) %%% check if centrosome is on right boundary
        cL = interp1(x(IL(1)-1:IL(end)+1),conc(IL(1)-1:IL(end)+1),xwrap(Left_Ind)); 
    else
        cL = zeros(1,size(Left_Ind,2)); 
        vec1 = xwrap(Left_Ind); 
        vec2 = vec1(vec1<x(end)-dx); 
        cL(1:size(vec2,2)) = interp1(x(IL(1):IL(end)),conc(IL(1):IL(end)),vec1(1:size(vec2,2)));
    end 
    
    
    %%% rest of cL(size(vec2,2)+1:end) remains 0. Since the density of MT ends
    %%% within a distance of dx from centrosome is 0, this has no consequence
    %%% on force calculation. This is correct if forces are assumed to be
    %%% coming from MT ends. If force is generated all along the MT length then
    %%% this will lead to error over a length scale of dx. This may be
    %%% significant since near the censtrosome ALL Mts will contribute to force
    %%% generation.

    fneg=dx_MT*trapz(cL(1:Ncen-1).*Pdfnew(1:Ncen-1));  %%%% force from negative branch

    %%%%% Right region split into two
    region1 = x(xcen_ind:end);    
    region2 = x(1:ceil(lenx/dx)+2); 
 
    IR = zeros(1,size(region1,2)+size(region2,2));
    IR(1:size(region1,2)) = xcen_ind:size(x,2); 
    IR(size(region1,2)+1:end) = 1:ceil(lenx/dx)+2; 

    %%%%% Ncen+size(Ip,2) =  centrosome coordinate in grid 2; 
    Right_Ind1 = size(Ip,2)+Ncen+1:size(xwrap,2);
    Right_Ind2 = 1:size(Ip,2) ;

    cR1= interp1(region1,conc(xcen_ind:end),xwrap(Right_Ind1));
    cR2= interp1(region2,conc(1:ceil(lenx/dx)+2),xwrap(Right_Ind2));

    fpos = dx_MT*trapz(cR1.*Pdfwrap(Right_Ind1)) +...
           dx_MT*trapz(cR2.*Pdfwrap(Right_Ind2)); 
    
 
%%%%%% CASE 5 : left branch in pervious periodic image.
elseif  (xlo<x(1) && xhi<x(end))  
    cs = 'case 5';                             

    In = find(xnew < x(1)); 
    xwrap = zeros(1,size(xnew,2));
    Pdfwrap = zeros(1,size(xnew,2)); 

    xwrap(1:size(xwrap,2)-size(In,2)) = xnew(In(end)+1:end);
    xwrap(size(xwrap,2)-size(In,2)+1:end) = x(end)+xnew(1:In(end));

    Pdfwrap(1:size(xwrap,2)-size(In,2)) = Pdfnew(In(end)+1:end);
    Pdfwrap(size(xwrap,2)-size(In,2)+1:end) = Pdfnew(1:In(end));

    %%%%%%% force calculation from right branch remains same as case 1 except 
    %%%%% that it gets shifted. 

    Ihi_grid1 = ceil(xhi/dx)+1; 
    IR = xcen_ind:Ihi_grid1+1; 
 
    %%%% Right branch ends at xhi which is the size(xwrap,2)-size(In,2)th point in
    %%%% xwrap  array. Centrosome index is Ncen-size(In,2)+1. We define branch
    %%%% from next point. 

    Right_Ind = (Ncen-size(In,2)+2):size(xwrap,2)-size(In,2); 
    cR = interp1(x(IR),conc(IR),xwrap(Right_Ind)); 
    fpos=dx_MT*trapz(cR.*Pdfwrap(Right_Ind));  %%%% force from negative branch


    %%%%% Left region split into two
    region1 = x(1:xcen_ind+2);  
    b_ind= floor(abs(xlo-2*dx)/dx);       %%% margin is built in here
    region2 = x(size(x,2)-b_ind-1:end);   %%% defined to make sure that the last point is x(end)

    IL = zeros(1,size(region1,2)+size(region2,2));
    IL(1:size(region1,2)) = 1:xcen_ind+2; 
    IL(size(region1,2)+1:end) = size(x,2)-size(region2,2)+1:size(x,2); 
 
    %%%%% Ncen-size(In,2) +1 =  centrosome coordinate in grid 2; 

    Left_Ind1 = 1:Ncen-size(In,2);
    Left_Ind2 = size(xwrap,2)-size(In,2)+1:size(xwrap,2);

    cL1= interp1(region1,conc(1:xcen_ind+2),xwrap(Left_Ind1));
    cL2= interp1(region2,conc(size(x,2)-size(region2,2)+1:size(x,2)),xwrap(Left_Ind2));

    fneg = dx_MT*trapz(cL1.*Pdfwrap(Left_Ind1)) +...
          dx_MT*trapz(cL2.*Pdfwrap(Left_Ind2)); 

end

rhs=fpos-fneg;

% end
    