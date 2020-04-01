%% Analysis of CuNi Sputtered Thin Film with RBS
% Steven Bopp, 2/8/2016

%% Integration for Compositional Information
Channel=linspace(500,540,41);
Counts=[4 3 7 3 9 6 11 26 34 52 111 183 342 465 478 439 338 298 237 219 196 207 231 278 399 444 442 382 313 270 218 161 86 46 20 8 6 11 4 2 2];

CuArea=5*trapz(linspace(500,520,41),Counts) % Units of Counts/keV
NiArea=5*trapz(linspace(520,540,41),Counts) % Multiply by five for 5 keV per channel


%% CuNi Double Peak Plot
% 5 keV per channel

plot(Channel,Counts,'p-k')
xlabel('Channel'); ylabel('Counts'); 
title('RBS Analysis of a CuNi Thin Film From Channel 500 to 540');
%% RBS Settings

%Group 1 CuNi  Tv=(1230) q=100k
%Spectrum:	RBS
%wilkens
%2/1/2016
%Geometry:	cornell
%Theta:	8
%Phi:	10
%Psi:	10
%Omega:	3
%BeamParameters:
%Energy:	4
%Z:	2
%Mass:	4
%Q:	2
%Charge:	1
%Current:	1
%MCAParameters:
%kev/ch:	5
%kevch0:	25
%StartingCh:	0
%FWHM:	25
