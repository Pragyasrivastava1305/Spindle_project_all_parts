function f = ellipse(b,ecc,N,sp_ang,lsp,Lm)
ms =3; lw =2;
a = b/sqrt(1-ecc^2);
th = linspace(0,2*pi,N);

r=@(th) b./sqrt(1-(ecc^2)*cos(th).^2); 
xunit = r(th).* cos(th) ;
yunit = r(th).* sin(th) ;
% for ith =1:size(th,2)
%     s(ith) = integral(r,0,th(ith)); 
% end

%%%%% angles of the two centrosomes.
%%%%% define the centrosome which is in 1st/4th quad to be cen1 and the one
%%%%% in 2nd/3rd quad to be cen2. 

if (cos(sp_ang)>0 || sp_ang == pi/2)
    phi1 = sp_ang; phi2 = mod(sp_ang + pi,2*pi); 
elseif (cos(sp_ang)<0 || sp_ang == 3*pi/2)
    phi1 = mod(sp_ang + pi,2*pi); phi2 = mod(sp_ang,2*pi);
end

%%% coordinates of centrosomes
C1x = lsp*cos(phi1); C1y = lsp*sin(phi1);
C2x = lsp*cos(phi2); C2y = lsp*sin(phi2);


    
f = figure(1);
plot(xunit,yunit,'ro','Markersize',ms);
    hold on; plot([C1x,C2x],[C1y,C2y],'sq-','Color','[1 0.6 0]','linewidth',lw,'Markersize',4*ms...
         ,'Markeredgecolor','k','Markerfacecolor','[1 0.6 0]');  
    hold on; 
    circle(C1x,C1y,Lm,lw,'b');   %%% ran discs 
    circle(C2x,C2y,Lm,lw,'b');  
    xlabel('x(\mum)','Fontsize',16)
    ylabel('y(\mum)','Fontsize',16)
    axis equal 
    xlim([-a-2  a+2]); ylim([-b-2,b+2])
    
    drawnow
end
