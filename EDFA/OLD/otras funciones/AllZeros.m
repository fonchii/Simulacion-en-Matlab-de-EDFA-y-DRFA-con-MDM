function z=AllZeros(f,xmin,xmax,N)
% Inputs :
% f : function of one variable
% [xmin - xmax] : range where f is continuous containing zeros
% N : control of the minimum distance (xmax-xmin)/N between two zeros
if (nargin<4)
    N=100;
end
dx=(xmax-xmin)/N;
x2=xmin;
y2=f(x2);
z=[];
options = optimoptions ('fsolve', 'Display', 'none');
for i=1:N
    x1=x2;
    y1=y2;
    x2=xmin+i*dx;
    y2=f(x2);
    if (y1*y2<=0)                              % Rolle's theorem : one zeros (or more) present
        z=[z,fsolve(f,(x2*y1-x1*y2)/(y1-y2), options)]; % Linear approximation to guess the initial value in the [x1,x2] range.
    end
end