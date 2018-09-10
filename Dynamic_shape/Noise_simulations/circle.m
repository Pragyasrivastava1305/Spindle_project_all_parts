function h = circle(x,y,r,m,col)
hold on
th = 0:pi/50:2*pi;
xunit = r .* cos(th) + x;
yunit = r .* sin(th) + y;
h = plot(xunit, yunit,'linewidth',m,'Color',col);
hold off