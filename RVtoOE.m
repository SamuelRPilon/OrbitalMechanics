clc
clear
close all

r_vect=[7500;0;0];
v_vect=[1;6;5];


r=norm(r_vect);
v=norm(v_vect);

r_hat=r_vect/r;

h_vect=cross(r_vect,v_vect);
h=norm(h_vect);

I=[1;0;0];
J=[0;1;0];
K=[0;0;1];

n=cross(K,h_vect)/norm(cross(K,h_vect));

mu=398600; %km^3/s^2

%%
%size and shape of the orbit
%semimajor axis: a
specificE=(v^2/2)-(mu/r);
a=-mu/(2*specificE);

%eccentricity: e
e_vect=((cross(v_vect,h_vect))/mu)-r_hat;
e=norm(e_vect);
e_hat=e_vect/e;

% energy of orbit 
E= v^2/2-mu/r;

%orbit Period
T= 2*pi*sqrt(a^3/mu);



%%
%display data
T = table([a],[e],[E],[T],'VariableNames',{'semimajor axis','eccentricity','Energy','Period'},'RowName',{});
disp(T)
