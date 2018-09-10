clear all
L=36; koffin = 0.1; koffout =0.02; dx = 0.01; kon =1;
x=-L:dx:L-dx;
Dvec = [0.001 0.01 0.1 1]; 
lvec = 5:0.1:15;

steady_state_sol = zeros(1,size(x,2));

for jD = 1:size(Dvec,2)
    jD
for jl = 1:size(lvec,2)
        
    D = Dvec(jD);  
    l = lvec(jl); 
    ld1 = sqrt(D/koffout); 
    ld2 = sqrt(D/koffin); 

    Num_fac = (ld1-ld2)*exp(2*l*(ld1+ld2)/ld1/ld2) - (ld1-ld2)*exp(2*L/ld1) - ...
            - (ld1+ld2)*exp(2*l/ld1) + (ld1+ld2)*exp((2*L/ld1) + 2*l/ld2); 

    A1 = exp((l+2*L)/ld1)*(exp(2*l/ld2)-1)*(koffout-koffin)*kon*ld1/koffout/koffin/Num_fac; 
    B1 = exp(l/ld1)*(exp(2*l/ld2)-1)*(koffout-koffin)*kon*ld1/koffout/koffin/Num_fac; 
    A2 = exp(l/ld2)*(exp(2*l/ld1)-exp(2*L/ld1))*(koffout-koffin)*kon*ld2/koffout/koffin/Num_fac;
    B2 = exp(l/ld2)*(exp(2*l/ld1)-exp(2*L/ld1))*(koffout-koffin)*kon*ld2/koffout/koffin/Num_fac;
    A3 = B1;
    B3 = A1;

    fun_reg1 = A1*exp(x/ld1) + B1*exp(-x/ld1) + kon/koffout; 
    fun_reg2 = A2*exp(x/ld2) + B2*exp(-x/ld2) + kon/koffin; 
    fun_reg3 = A3*exp(x/ld1) + B3*exp(-x/ld1) + kon/koffout; 

    %%%%% Define piecewise function
    steady_state_sol(x <= -l) = fun_reg1(x <= -l); 
    steady_state_sol(x > -l & x <= l) = fun_reg2(x > -l & x <= l); 
    steady_state_sol(x > l) = fun_reg3(x > l); 
    I = find(steady_state_sol < 0.9*steady_state_sol(end)); 

    lsim(jD,jl) = 0.5*size(I,2)*dx; 
    
    figure(1)
    if jl == 61
        plot(x,steady_state_sol,'linewidth',1.5); drawnow; hold on
    end 

end 
figure(2)
plot(lvec,lsim(jD,:),'linewidth',1.5);  drawnow;  hold on 

end 
figure(1)
plot(x,0.9*steady_state_sol(end)*ones(1,size(x,2)),'-.k','linewidth',1.5)
grid on; 
xlim([0 L])

figure(2)
plot(lvec,11*ones(1,size(lvec,2)),'-.k','linewidth',1.5)
grid on 










