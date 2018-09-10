
Dvec = [0.001 0.01 0.1 1];
lvec = 2:0.5:10; 

% for jD = 1:size(Dvec,2)-1
jD = 1
for jl = 1:size(lvec,2)
    l = lvec(jl); 
    D = Dvec(jD); 
    load(['D=',num2str(D),'l=',num2str(l),'.mat'])
    
    
    plot(x,cL_arr(:,end));hold on
    I = find(cL_arr(:,end) <= 0.9*cL_arr(1,end)); 
    Lout(jD,jl) = 0.5*dx*size(I,1);

end
% end