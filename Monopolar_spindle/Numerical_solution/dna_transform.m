function [xdna,cdna,h] = dna_transform(x,xN,conc,dx,koffin)
%%%  define x coordinate and LGN concentration in dna frame

N =  size(x,2); 
xnew = x-xN; 
 
%%%% bring 0 in middle
ind_vec = find(xnew >= 0.5*x(end));

xdna = linspace(-0.5*x(end)+dx,0.5*x(end),N);
cdna = koffin*circshift(conc,size(ind_vec,2)); 

h = plot(xdna,cdna,'linewidth',1.5);
hold on 
scatter(0,0)

Mat = [xdna' cdna];
dlmwrite(['steady_LGN_D=0.01_vmean=10.22_eps=0.csv'],Mat,'delimiter','\t')

