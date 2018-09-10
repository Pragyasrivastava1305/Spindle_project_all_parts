Df = 1;
L = 72;
dx = 0.1;
x = 0:dx:L;
N = size(x,2);

delt = 0.001; tmax = 1;
time = 0:delt:tmax;
Niter =  size(time,2)
 
cL_old = exp(-(x-36).^2/2);

for iter =1:Niter-1
    
    der2 = (circshift(cL_old,-1)+circshift(cL_old,1)-2*cL_old)/dx/dx;
    cL_new = cL_old+delt*Df*der2;
   
    cL_old = cL_new;
    
    plot(x,cL_new,'linewidth',2,'color','r')
    drawnow
end