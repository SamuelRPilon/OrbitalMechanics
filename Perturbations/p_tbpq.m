function [sdo]=p_tbpq(t,si)


mu = 3.986e5; %km^3 s^-2
re1= 6378;
j2= 1.08263e-3;
r=norm(si(1:3)); 
ac= -(3*j2*mu*(re1^2))/(2*r^4);

sdo(1)= si(4);
sdo(2)= si(5);
sdo(3)= si(6);

sdo(4)= -mu*si(1)/r^3 + ac*(1-5*(si(3)/r)^2) * (si(1)/r);
sdo(5)= -mu*si(2)/r^3 + ac*(1-5*(si(3)/r)^2) * (si(2)/r);
sdo(6)= -mu*si(3)/r^3 + ac*(1-5*(si(3)/r)^2) * (si(3)/r);

sdo = sdo'; 
end