function [dist,Pts] = closest_dist(r, th, p1,p2,varargin)
%%%%% INPUT:  shape r of the ellipse, th:  angular grid,
%%%%% p1 and p2 contain x-y coordinate of the two end points of metaphase
%%%%% plate. 
%%%% calculates closest distance from the points on periphery to the line
%%%% segment defined by the two end points p1 and p2. 
%%%% ALGORITHM :  for a given point, drop a perpendicular to the line
%%%% segment, if the end point of this perpendicular is on the line
%%%% segment, assign the length to dist. If it falls out of the line
%%%% segment, the assign minimum of the distances to two end point to dist.

%%% for troubleshooting and to run it as a code uncomment the following
% a = 20;  b =10; lA=9;  ang1 = pi/4; 
% shape = @(th) a./sqrt(a^2/b^2 - (a^2/b^2-1)*cos(th).^2);
% N1x = lA*cos(ang1); N1y = lA*sin(ang1);
% N2x = lA*cos(ang1+pi); N2y = lA*sin(ang1+pi);
%p1 = [N1x N1y]; p2 = [N2x N2y]; 
%r =shape; 

flag = ismember('flag',varargin); 

%%%% find the coordinates of intersection point 
x1 = p1(1); x2= p2(1); y1 = p1(2); y2= p2(2); 
len = sqrt((x1-x2)^2 + (y1 - y2)^2); 
xp = r(th).*cos(th); yp = r(th).*sin(th);
xR = (xp.*x1^2 + x2*(x2.*xp + (y1-y2)*(y1-yp)) + x1*(-2*x2*xp-(y1-y2).*(y2-yp)))/len/len;
yR = ((x1-x2)*(-x2*y1 + xp*(y1-y2) + x1*y2) + yp*(y1-y2)^2)/len/len; 
%Pt = [xR; yR]';   %%%% N times 2 matrix where N is size of th matrix

%%% distance between ends of metaphase plate and the point of perp drop
%%% both of these should be less than len for the perp drop on the segment.
N1_2_R = sqrt((x1-xR).^2 + (y1-yR).^2); 
N2_2_R = sqrt((x2-xR).^2 + (y2-yR).^2); 

%%% distance between the point on periphery and the two ends of the plate
p_to_N1 = (x1-xp).^2 + (y1-yp).^2;    
p_to_N2 = (x2-xp).^2 + (y2-yp).^2;
min_dist_from_points = sqrt(min(p_to_N1,p_to_N2)); 

dist = min_dist_from_points; 
%plot(dist);  hold on 
%%%% set dist to the length of perpendicular if the point of contact drops
%%%% between the end points. 
perp_len = sqrt(((xp.*(y2-y1) + yp.*(x1-x2)+x2*y1 - x1*y2).^2)/len/len); 
Ivec =  find((N1_2_R <= len) & (N2_2_R <= len)); 

dist(Ivec) = perp_len(Ivec); 
%plot(dist)

%%%% get the point of closest distance :  first set it one of the end
%%%% points based on which one is closest, and then for the perp drop on
%%%% the line segment, set it to xR,yR. 

Pts(p_to_N1 <= p_to_N2,1) = x1; 
Pts(p_to_N1 <= p_to_N2,2) = y1; 
Pts(p_to_N1 > p_to_N2,1) = x2; 
Pts(p_to_N1 > p_to_N2,2) = y2;

Pts(Ivec, 1) = xR(Ivec); 
Pts(Ivec, 2) = yR(Ivec); 


% %%%%% plot if flag is true
if flag == true
%  figure; 
    plot(xp,yp,'ro','Markersize',3); hold on;    
    for ith =1:size(th,2)
    plot([xp(ith),Pts(ith,1)],[yp(ith),Pts(ith,2)],'color',[0.2 0.7,0])    
    end
    plot([x1,x2],[y1,y2],'linewidth',1.5); 
    drawnow
    axis equal
    hold off
 end


end

