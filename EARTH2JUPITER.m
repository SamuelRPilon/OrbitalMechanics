
%% SAMUEL R PILON
% AE 313 - EXAM #2 
% For personal use of SAMUEL PILON only. Not to be distributed

clc; 
close; 
clear; 
%% QUESTION STATEMENT 
%% Use A Hohmann Transfer To Get from Earth To Juipter

%% GIVEN INFORMATION
%
EarthParkOrbit = 285; %km       %Parking Orbit around Earth 
a = 778.6*10^6; %km             %Semi-Major Axis of a Heliocentric Hohmann Transfer Orbit
JupiterParkOrbit = 80000; %km   %Parking Capture Orbit around Jupiter
RadiusJ = 71490; %km            %Equatorial Radius of Jupiter 
RadiusE = 6378;  %km            %Equatorial Radius of Earth
muJ = 126686000; %km3s-2        %Standard gravitational Parameter of Jupiter
muS = 132.7*10^9; %km3s-2       %Standard gravitational Parameter of the Sun
muE = 398600; %km3s-2           %Standard gravitational Parameter of Earth
ESD = 147.6*10^6; %km           %Distance from the Earth To the Sun 
JSD = 778600000;  %km           %Distance from Jupiter To the Sun
RE = EarthParkOrbit + RadiusE;  %Total Distance from Earth to Parking orbit
RJ = JupiterParkOrbit + RadiusJ;%Total Distance from Jupiter to Parking orbit


%% 1. What is the ∆v required to leave Earth’s orbit for Jupiter?

% From Earth To Jupiter 
Ve = sqrt(muS/ESD);                          %Velicity of Earth
EnergyT = -muS/(ESD+JSD);                    %Energy of The Total Transfer
V1 = sqrt(2 * ((muS / ESD) + EnergyT));      %Departure Velocity of Earth
V_infe = abs(V1 - Ve);                       %Transfer Velocity Leaving Earth to Jupiter

% Departure from Earth 
EnergyinfE = (V_infe^2)/2;                   %Specific Energy of Hyperbolic Trajectory
Vhype = sqrt(2 * ((muE/RE) + EnergyinfE));   %Hyperbolic Velocity Injection
Vpark = sqrt(muE/RE);                        %Velocity of The Parking Orbit 
delVBoost = abs(Vhype - Vpark);              %Delta V of The Boost Needed For The Target 

fprintf('The required ∆v to leave Earth’s orbit for Jupiter is = %.3f\n',delVBoost);
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('                  These Are The Values Used To Calculate This Question');
fprintf('\n');
fprintf('\n');
p=table(Ve, EnergyT, V1,V_infe,EnergyinfE,Vhype,Vpark,delVBoost);
disp(p);
fprintf('-----------------------------------------------------------------------------------------------\n');

%% 2. What is the departure hyperbola eccentricity and what is the departure angle, η ?

%The radius to periapse of the departure hyperbola is the radius of 
%the earth plus the altitude of the parking orbit,

e = 1 + (RE * V_infe^2)/muE;                 %Eccentricity of the Departure Hyperbola
n = acosd(-1/e);                             %Departure Angle η


fprintf('The Departure Hyperbola Eccentricity is %.3f and The Departure Angle is %.3f°\n',e,n);
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('                  These Are The Values Used To Calculate This Question');
fprintf('\n');
fprintf('\n');
o=table(e,n,V_infe);
disp(o);
fprintf('-----------------------------------------------------------------------------------------------\n');

%% 3. If the spacecraft stops at Jupiter, as planned, what is the total mission ∆v,
%     from Earth parking orbit through planetary capture at Jupiter?

% Arivial at Jupiter
V2 = sqrt(2*((muS/JSD)+EnergyT));            %Velocity of the Arrivial at the Target
Vj = sqrt(muS/JSD);                          %Velocty of Jupiter
VinfT= abs(Vj - V2);                         %Injection Velocty at Jupiter
EinfT = (VinfT^2)/2;                         %Energy of the Transfer
VHypeT = sqrt(2*((muJ/RJ)+EinfT));           %Hyperbolic Velocity Injection
VparkT = sqrt(muJ/RJ);                       %Velocity of the Parking Orbit 
delVRetro = abs(VparkT - VHypeT);            %Delta V of the Reverese Thrust Needed to Slow Down
    
DeltaVmission = delVBoost + delVRetro; 
fprintf('The  ∆v of the Total Mission  is = %.3f\n',DeltaVmission);
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('                  These Are The Values Used To Calculate This Question');
fprintf('\n');
fprintf('\n');
j=table(V2, Vj,VinfT,EinfT,VHypeT,VparkT,delVRetro);
disp(j);
fprintf('-----------------------------------------------------------------------------------------------\n');

%% 4. If the spacecraft malfunctions and does not burn, what is the ∆v about the Sun
%     considering the gravitational interaction with Jupiter?
    
a1 = .5*(ESD+JSD);                            %Semimajor Axis of the Hohmann Transfer Ellipse
V1v = sqrt(muS*((2/JSD)-(1/a1)));             %Spacecraft Velocity Upon Arivial in Jupiters Sphere of Influence
Vinf4 = abs(Vj-V1v);                          %Speed of the Escape Hyperbola 
e1 = 1 + (RJ*Vinf4^2)/muJ;                    %Eccentricity of the Hyperbolic Swing 
del = 2*asind(1/e1);                          %Turn angle 

DeltaVJupiterG = 2*Vinf4*sind(del/2); 
fprintf('The  ∆v of Jupiters Gravitational Interaction is = %.3f\n',DeltaVJupiterG);
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('                  These Are The Values Used To Calculate This Question');
fprintf('\n');
fprintf('\n');
q=table(a1,V1v,Vinf4,e1,del);
disp(q);
fprintf('-----------------------------------------------------------------------------------------------\n');

%% 5. Plot the heliocentric ∆v from a planetary flyby for a range of aiming radii, ∆, from
%     100,000 km to 2,000,000 km (in steps of at most 10,000 km).

delta = 100000:10000:2000000;               %Spans the aiming radii from 100000 to 2000000
Q = (2*muJ)/(Vinf4^2);                      %Calculates part for ease of format
rp = -(Q + sqrt(Q^2 + 4*delta.^2))/2;       %Radius of perigee with span Delta
Einf5 = Vinf4^2/2;                          %Energy of the span per Delta
h = rp.*sqrt(Vinf4^2+((2*muJ)./rp));        %Angular momentum per span delta
e5 = sqrt(1+((2*Einf5.*h.^2)/muJ^2));       %Eccentricity of departure orbit 
deltaV5 = (2*Vinf4)./e5;                    %Delta V

plot(delta,deltaV5,'r--o','LineWidth',.1);
title('Heliocentric ∆v From a Planetary Flyby for a Range of Aiming Radii ∆'); 
xlabel('Aiming Radii "∆"');
ylabel('\Delta V')
grid on 

fprintf('\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('                  These Are The Values Used To Calculate This Question');
fprintf('\n');
fprintf('\n');
w = table(Q,Einf5);
disp(w);
fprintf('-----------------------------------------------------------------------------------------------\n');

