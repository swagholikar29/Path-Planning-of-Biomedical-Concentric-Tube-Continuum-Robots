%% Image Analyzer Script
% goes through series of images and applies the filter to get the tube
% then uses circle fit on the tube points

th = linspace(0,2*pi,20)';

img = imread('cropped.jpg');
[bw, rgb] = maskBlackTubes(img);

[y, x] = find(bw);
[xc,yc,Re,a] = circfit(x,y);
xe = Re*cos(th)+xc; ye = Re*sin(th)+yc;

plot(x,y,'o',[xe;xe(1)],[ye;ye(1)],'-.',R*cos(th),R*sin(th));
title(' measured fitted and true circles')
legend('measured','fitted','true')
text(xc-R*0.9,yc,sprintf('center (%g , %g );  R=%g',xc,yc,Re))
xlabel x, ylabel y 
axis equal