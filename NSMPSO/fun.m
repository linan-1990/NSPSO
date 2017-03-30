function [f1,f2]=fun(x)
f1=sum((x).^2);
f2=sum((x-2).^2);

%function [f1,f2]=fun(x1,x2,x3)
%f1=-10*exp(-0.2*sqrt(x1^2+x2^2))-10*exp(-0.2*sqrt(x2^2+x3^2));
%f2=(abs(x1))^0.8+5*sin(x1^3)+(abs(x2))^0.8+5*sin(x2^3)+(abs(x3))^0.8+5*sin(x3^3);