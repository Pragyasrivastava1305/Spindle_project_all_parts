function [rhs] = Rates(kin,kout,x,dx,xN,l,L)
%%%%% Given the position of nucleus and range l this function calculates rates 
%%%%% If the nucleus position falls in between two grid
%%%%% points, then a fractional contribution is added at the grid points
%%%%% next to edges x1 and x2. l < L/2 is required. 

%%% Obtain left and right edges of the inhibited range
x1=xN-l; x2=xN+l; 

if (x1>=0 && x2<=L)
    x1p=x1; x2p=x2;
    
elseif x1<0
    x1p=L+x1; x2p=x2;

elseif x2>L
   x1p=x1; x2p=x2-L;
end


Lind=ceil(x1p/dx)+1; Rind=floor(x2p/dx)+1; 
%%% indices of left and right edges. 1 added because first point is x=0.

f=mod(x2p,dx)/dx;     
%%% the fraction of length that may be between two grid points on chemical
%%% grid

rhs=zeros(1,size(x,2)); 

%%%%% Take into account periodicity
if (x1>=0 && x2<=L)        %%%% When inhibited range is inside the chemical grid

rhs(Lind:Rind)=kin;        %%%% points inside the DNA region

rhs(1:Lind-2)=kout;        %%%% points outside the DNA region
rhs(Rind+2:end)=kout;

rhs(Lind-1)=kin+(1-f)*(kout-kin);   %%%%% linear interpolation just
rhs(Rind+1)=kin+f*(kout-kin);       %%%%% outside the edge
    
    
elseif (x1<0 || x2>L)       %%%% When left edge < 0  or right edge > L
    
rhs(1:Rind)=kin;  rhs(Lind:end)=kin;

rhs(Rind+2:Lind-2)=kout; 

rhs(Lind-1)=kin+(1-f)*(kout-kin);

rhs(Rind+1)=kin+f*(kout-kin);




end






end

