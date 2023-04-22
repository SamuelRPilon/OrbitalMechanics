function [sdo]=p_tbp(t,si)


mu = 3.986e5; %km^3 s^-2
r=norm(si(1:3)); 
sdo(1)= si(4);
sdo(2)= si(5);
sdo(3)= si(6);
sdo(4)= -mu*si(1)/r^3;
sdo(5)= -mu*si(2)/r^3;
sdo(6)= -mu*si(3)/r^3;

sdo = sdo'; 
end