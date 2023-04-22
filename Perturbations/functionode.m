clc 
clear all

r = 6378; %km
G = 6.674e-20; %km^3kg^-1s^-2
mu = 3.986e5; %km^3 s^-2
J2 = 1.08263e-3;
tspan =[0,9000];
r1 = [7500 0 0];
v1 = [1 6 5];
y=[7500 0 0 1 6 5];
k = [0 0 1];
h1 = cross(r1,v1);
i1 = dot(h1,k);

[tout,yout]=ode45(@(tout, yout) p_tbp(tout, yout),tspan,y)
[x y z]= sphere;

x=x*r; 
y=y*r;
z=z*r;
surf(x,y,z)
hold on 
title('Two Body Orbit')

plot3(yout(:,1), yout(:,2), yout(:,3),'color', 'r')
hold off 
axis equal
