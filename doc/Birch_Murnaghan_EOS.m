%% VASP Homework Assignment Due 11-20-2015
clear all, clc
%% Birch Murnaghan eqn. of State
% Re-fit with cftool and the Birch Murnaghan equation of state:
% E+(9*V*B/16)* (((A/x)^2-1)^3 * D+((A/x)^2 -1)^2*(6-4*(A/x)^2 ))
%% Diamond FCC
E0= [-.19639643E+02, -.20138523E+02, -.20331547E+02, -.20277141E+02, -.20027958E+02];
x=[3.3:0.1:3.7]; xx=[3.3:.001:3.7];
yy=spline(x,E0,xx);        % Create spline interpolate for the lattice constant

figure(1);
plot(x,E0,'m*',xx,yy,'g'); % Plot spline interpolate and Energy with Lattice Parameter
xlabel('Lattice Parameter in Angstroms'); ylabel('Energy in eV');
title('Energy Minimization for Diamond Cubic Carbon');
% Pentagram p; hexagram h; diamond d; square s;

indexmin=find(min(yy) == yy);               % Define indexmin
xmin = xx(indexmin); ymin = yy(indexmin);   % Calculate minimum values
A0 = xmin; A0error = (3.526-A0)/3.526 *100; % Percent error
fprintf('Lattice Constant for diamond:%g Angstroms \n',A0) % Display lattice constant and error from expected value
fprintf('Percent error between calculated and experimental lattice parameters for diamond:%g \n',A0error)

%% Copper FCC
figure(3); hold on;
E0=[-.31504565E+01 -.34863921E+01 -.36698485E+01 -.37354998E+01 -.37362657E+01 -.37072386E+01 -.36671771E+01];
x=[3.3 3.4 3.5 3.6 3.65 3.75 3.8]; xx=[3.3:.001:3.8];
yy=spline(x,E0,xx); % Create spline interpolate for the lattice constant
plot(x,E0,'m*',xx,yy,'g'); title('Copper FCC Energy Minimization as a Function of Lattice Parameter');
xlabel('Lattice Parameter'); ylabel('Energy in eV');

indexmin=find(min(yy) == yy);               % Define indexmin
xmin = xx(indexmin); ymin = yy(indexmin);   % Calculate minimum values
A0 = xmin; A0error = (3.615-A0)/3.615 *100; % Percent error
fprintf('Lattice Constant for FCC Cu:%g Angstroms \n',A0) % Display lattice constant and error from expected value
fprintf('Percent error between calculated and experimental lattice parameters for FCC Cu:%g \n',A0error)

%% Copper SC
%With Avogadro I found that the lattice parameter should be about 2.604 angstroms.

%figure(4);
E0=[-.31037012E+01 -.32108047E+01 -.31709793E+01 -.30421292E+01 -.28702942E+01 -.26680372E+01 -.16849119E+01 -.15256076E+01 -.13754908E+01 -.12413620E+01 -.11193245E+01 -.10086932E+01];
x=[2.3 2.4 2.5 2.6 2.7 2.8 3.3 3.4 3.5 3.6 3.7 3.8]; xx=[2.3:.001:3.8];
yy=spline(x,E0,xx); % Create spline interpolate for the lattice constant
plot(x,E0,'r*',xx,yy,'b'); title('Copper SC Energy Minimization as a Function of Lattice Parameter');
xlabel('Lattice Parameter'); ylabel('Energy in eV'); axis([ 2 4 -3.8 -2.5]);

indexmin=find(min(yy) == yy);               % Define indexmin
xmin = xx(indexmin); ymin = yy(indexmin);   % Calculate minimum values
A0 = xmin;
fprintf('Lattice Constant for Simple Cubic Cu:%g Angsroms which was expected to be about 2.604 Angstroms \n',A0) % Display lattice constant and error from expected value

%% Copper BCC e.q. lattice parameter should be about 2.9 angstroms

E0=[-.34304587E+01 -.36437209E+01 -.36898924E+01 -.36342168E+01 -.35102234E+01 -.33427409E+01 -.31517294E+01];
x=[2.7:0.1:3.3]; xx=[2.7:.001:3.3];
yy=spline(x,E0,xx); % Create spline interpolate for the lattice constant
plot(x,E0,'k*',xx,yy,'m'); title('Copper FCC, BCC and SC Energy Minimization as a Function of Lattice Parameter');
xlabel('Lattice Parameter'); ylabel('Energy in eV');

indexmin=find(min(yy) == yy);               % Define indexmin
xmin = xx(indexmin); ymin = yy(indexmin);   % Calculate minimum values
A0 = xmin;
fprintf('Lattice Constant for BCC Cu:%g Angstroms which was expected to be about 2.9 Angstroms \n',A0) % Display lattice constant and error from expected value
