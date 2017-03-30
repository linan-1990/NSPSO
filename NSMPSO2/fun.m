function [f1,f2]=fun(x1,x2)
f1=x1^2+x2^2;
f2=(x1-2)^2+(x2-2)^2;
% function [f1,f2]=fun(x,y)
% global Ti Vi PMV Energy;
% f1 = (interp2(Ti,Vi,PMV,x,y,'spline'))^2;
% f2 = interp2(Ti,Vi,Energy,x,y,'spline');