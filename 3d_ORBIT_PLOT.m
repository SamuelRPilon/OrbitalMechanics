%% Samuel R Pilon
% AE 313 
% For personal use of Samuel R Pilon only. Not to be distributed.

%% Starts the Function
function orbit
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%{
  This function computes the orbit of a spacecraft.
 
  It also plots the orbit and computes the times at which the maximum
  and minimum radii occur and the speeds at those times.
 
  %}
% -------------------------------------------------------------------
% Clear the command window, clear all variables, and close all figures.
clc; close all; clear all
 
hours = 3600;
G     = 6.6742e-20;
 
%...Input data:
%   Earth:
m1 = 5.974e24;
R  = 6378;
m2 = 1000;
 
r0 = [8000 0 6000];
v0 = [0 7 0];
 
t0 = 0;
tf = 4*hours;
%...End of the 
 
 
%...Numerical integration:
mu    = G*(m1 + m2);
y0    = [r0 v0]';
[t,y] = rkf45(@rates, [t0 tf], y0);
 
%...Output the results:
output
 
return
 
% ~~~~~~~~~~~~~~~~~~~~~~~~
function dydt = rates(t,f)
% ~~~~~~~~~~~~~~~~~~~~~~~~
% This function calculates the acceleration vector 
% ------------------------
x    = f(1);
y    = f(2);
z    = f(3);
vx   = f(4);
vy   = f(5);
vz   = f(6);
 
r    = norm([x y z]);
 
ax   = -mu*x/r^3;
ay   = -mu*y/r^3;
az   = -mu*z/r^3;
    
dydt = [vx vy vz ax ay az]';    
end %rates
 
 
% ~~~~~~~~~~~~~
function output
% ~~~~~~~~~~~~~
%This function computes the maximum and minimum radii
% -------------
for i = 1:length(t)
    r(i) = norm([y(i,1) y(i,2) y(i,3)]);
end
 
[rmax imax] = max(r);
[rmin imin] = min(r);
 
v_at_rmax   = norm([y(imax,4) y(imax,5) y(imax,6)]);
v_at_rmin   = norm([y(imin,4) y(imin,5) y(imin,6)]);
 
%...Output to the command window:
fprintf('\n\n--------------------------------------------------------\n')
fprintf('\n Earth Orbit\n')
fprintf(' %s\n', datestr(now))
fprintf('\n The initial position is [%g, %g, %g] (km).',...
                                                     r0(1), r0(2), r0(3))
fprintf('\n   Magnitude = %g km\n', norm(r0))
fprintf('\n The initial velocity is [%g, %g, %g] (km/s).',...
                                                     v0(1), v0(2), v0(3))
fprintf('\n   Magnitude = %g km/s\n', norm(v0))
fprintf('\n Initial time = %g h.\n Final time   = %g h.\n',0,tf/hours) 
fprintf('\n The minimum altitude is %g km at time = %g h.',...
            rmin-R, t(imin)/hours)
fprintf('\n The speed at that point is %g km/s.\n', v_at_rmin)
fprintf('\n The maximum altitude is %g km at time = %g h.',...
            rmax-R, t(imax)/hours)
fprintf('\n The speed at that point is %g km/s\n', v_at_rmax)
fprintf('\n--------------------------------------------------------\n\n')
 
%...Plot the results:
%   Draw the planet
[xx, yy, zz] = sphere(100);
surf(R*xx, R*yy, R*zz)
colormap(light_gray)
caxis([-R/100 R/100])
shading interp
 
%   Draw and label the X, Y and Z axes
line([0 2*R],   [0 0],   [0 0]); text(2*R,   0,   0, 'X')
line(  [0 0], [0 2*R],   [0 0]); text(  0, 2*R,   0, 'Y')
line(  [0 0],   [0 0], [0 2*R]); text(  0,   0, 2*R, 'Z')
 
%   Plot the orbit, draw a radial to the starting point
%   and label the starting point (o) and the final point (f)
hold on
plot3(  y(:,1),    y(:,2),    y(:,3),'k')
line([0 r0(1)], [0 r0(2)], [0 r0(3)])
text(   y(1,1),    y(1,2),    y(1,3), 'o')
text( y(end,1),  y(end,2),  y(end,3), 'f')
 
%   Select a view direction (a vector directed outward from the origin) 
view([1,1,.4])
 
%   Specify some properties of the graph
grid on
axis equal
xlabel('km')
ylabel('km')
zlabel('km')
 
% ~~~~~~~~~~~~~~~~~~~~~~~
function map = light_gray
% ~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function creates a color map for displaying the planet as light
  gray with a black equator.
%}
% -----------------------
r = 0.8; g = r; b = r;
map = [r g b
       0 0 0
       r g b];
end %light_gray
 
end %output
 
end %orbit
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [tout, yout] = rkf45(ode_function, tspan, y0, tolerance)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{ 
  This function uses the Runge-Kutta-Fehlberg 4(5) 
%}
% ---------------------------------------------------------------
 
a = [0 1/4 3/8 12/13 1 1/2];
 
b = [    0          0          0          0         0
        1/4         0          0          0         0
        3/32       9/32        0          0         0
     1932/2197 -7200/2197  7296/2197      0         0
      439/216      -8      3680/513   -845/4104     0
       -8/27        2     -3544/2565  1859/4104  -11/40];
 
c4 = [25/216  0  1408/2565    2197/4104   -1/5    0  ];
c5 = [16/135  0  6656/12825  28561/56430  -9/50  2/55]; 
 
if nargin < 4
    tol  = 1.e-8;
else
    tol = tolerance;
end
 
t0   = tspan(1);
tf   = tspan(2);
t    = t0;
y    = y0;
tout = t;
yout = y';
h    = (tf - t0)/100; % Assumed initial time step.
 
while t < tf
    hmin = 16*eps(t);
    ti   = t;
    yi   = y;
    %...Evaluate the time derivative(s) at six points within the current
    %   interval:
    for i = 1:6
        t_inner = ti + a(i)*h;
        y_inner = yi;
        for j = 1:i-1
            y_inner = y_inner + h*b(i,j)*f(:,j);
        end
        f(:,i) = feval(ode_function, t_inner, y_inner);
    end
 
    %...Compute the maximum truncation error:
    te     = h*f*(c4' - c5'); % Difference between 4th and
                              % 5th order solutions
    te_max = max(abs(te));    
   
    %...Compute the allowable truncation error:
    ymax       = max(abs(y));
    te_allowed = tol*max(ymax,1.0);
    
    %...Compute the fractional change in step size:
    delta = (te_allowed/(te_max + eps))^(1/5);
     
    %...If the truncation error is in bounds, then update the solution:
    if te_max <= te_allowed
        h     = min(h, tf-t);
        t     = t + h;
        y     = yi + h*f*c5';      
        tout  = [tout;t];
        yout  = [yout;y'];
    end
    
    %...Update the time step:
    h  = min(delta*h, 4*h);
    if h < hmin
        fprintf(['\n\n Warning: Step size fell below its minimum\n'...
                 ' allowable value (%g) at time %g.\n\n'], hmin, t)
        return
    end  
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end