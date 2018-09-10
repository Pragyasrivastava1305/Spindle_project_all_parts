function [rhs]=rk4_RDM(c_old,lD,kon,koff,dt,dx)

dim=2;


Der2=(circshift(c_old,-1,dim)+circshift(c_old,1,dim)-2*c_old)/(dx^2);
sol1=koff.*((lD.^2).*Der2+(kon./koff)-c_old);
k1=dt*sol1; arg1=c_old+0.5*k1;

Der2=(circshift(arg1,-1,dim)+circshift(arg1,1,dim)-2*arg1)/(dx^2);
sol2=koff.*((lD.^2).*Der2+(kon./koff)-arg1);
k2=dt*sol2;   arg2=c_old+0.5*k2;

Der2=(circshift(arg2,-1,dim)+circshift(arg2,1,dim)-2*arg2)/(dx^2);
sol3=koff.*((lD.^2).*Der2+(kon./koff)-arg2);
k3=dt*sol3;  arg3=c_old+k3;


Der2=(circshift(arg3,-1,dim)+circshift(arg3,1,dim)-2*arg3)/(dx^2);
k4=dt*koff.*((lD.^2).*Der2+(kon./koff)-arg3); 

rhs=k1/6+k2/3+k3/3+k4/6; 


end

