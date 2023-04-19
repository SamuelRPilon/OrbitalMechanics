%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Samuel Pilon 
%Class: AE313 
%Date: 11/28/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
clear all 
close all 


PO = 320;%km
a = 227.9*10^6;%km
altitude_Mars = 800; %km
mu_mars = 42828; %km^3s^-2
mu_sun = 1.327*10^11; %km^3s^-2
mu_earth = 398600; %km^3s^-2
RMars = 227.9*10^6; %km 
REarth = 149.6*10^6; %km 
EQR_MARS = 3389.5; %km 
EQR_EARTH = 6378; %km 
ISP = 385; %seconds
go= 9.81*10^-3; 
RE= 6663;%EQR_EARTH + PO; 
RM = 3396;%EQR_MARS + altitude_Mars; 


%% From Earth To Target 
Ve = sqrt(mu_sun/REarth);           % velicity of Earth
EnergyT = -mu_sun/(REarth+RMars);     % energy of the totoal transfer
V1 = sqrt(2*((mu_sun/REarth)+EnergyT)); %departurn velocity of earth
vinfe = abs(V1-Ve);                 % transfer velocity leaving earth to target

Vtarget = sqrt(mu_sun/RMars);       %velocty of mars
V2 = sqrt(2*((mu_sun/RMars)+EnergyT));  %vlcity of the arrivial at target

%% Departure from Earth 
EnergyinfE = (vinfe^2)/2;   %specific energy of hyperbolic trajectory
vhype = sqrt(2*((mu_earth/RE)+EnergyinfE)); %hyperbolic velicty injection
Vpark = sqrt(mu_earth/RE);      %velocity of the parking orbit 

delVBoost = abs(vhype - Vpark); % delta v of the boost needed for the target 

%% Arivial at Mars 
Vm = sqrt(mu_sun/RMars);        %velocty of mars
VinfT= abs(Vm - V2);            %injection velocty at mars
EinfT = (VinfT^2)/2;            
VHypeT = sqrt(2*((mu_mars/RM)+EinfT)); %hyperbolic velicty injection
VparkT = sqrt(mu_mars/RM);              %velocty of the parking orbit 
delVRetro = abs(VparkT - VHypeT);       %delta v of the reverese thrust needed to slow dwn 


DeltaVmission = delVBoost + delVRetro; 
fprintf('delta v mission %f \n',DeltaVmission);

%% Determine the amout of propellant Required 

massper = 1 - exp(-DeltaVmission/(ISP*go));
massper_per = massper*100;

fprintf('The total mass percentage of the mission  is %f %',massper_per);
