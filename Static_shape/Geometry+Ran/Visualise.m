A= 2; 
ellipse(b,A,3600,pi/4,9,15); 
hold on 

[dist,Pt]= closest_dist(shape,th,[N1x N1y],[N2x,N2y]); 

xunit = shape(th).*cos(th); 
yunit = shape(th).*sin(th); 

for ith =1:size(dist,2)
    if dist(ith) < 12
    plot([xunit(ith),Pt(ith,1)],[yunit(ith),Pt(ith,2)],'color',[0.1 0.8,0])    
    hold on; drawnow; 
    end
    
end