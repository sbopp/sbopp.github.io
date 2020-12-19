%%
echo on
%-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
%-:-:-:-:-:-:-:-:-:-:-Symbolic_to_Vector_Piecewise_Function.m-:-:-:-:-:-:-:-:-:-:-
%-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
%-:-:-:-:-:-:-:-:-:Steven E. Bopp, Materials Science & Engineering:-:-:-:-:-:-:-:-
%-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-December 12, 2020-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
%-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-#PandemicMath-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
%-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
echo off

% Use this code to turn a symbolic function that defines a piecewise vector into a vector of type double
% I'm using it to draw the space for a finite quantum well

%%

clear all; close all; clc;

% This structure is defined as: barrier  |    well   |  barrier  |    well   |  barrier
%                               Region 1    Region 2    Region 3    Region 4    Region 5

r1V = 9; % Region 1 potential
r2V = 0; % Region 2 potential
r3V = 9; % Region 3 potential
r4V = 0; % Region 4 potential
r5V = 9; % Region 5 potential

r2w = 1; % Region 1 width
r2c = 0; % Region 1 center

syms x
y = piecewise(...
    x<-r2w/2, r1V, ... 
    x>-r2w/2 & x<r2w/2, r2V, ...
    x>r2w/2 & x<1+r2w/2, r3V, ...
    x>1+r2w/2 & x<2+r2w/2, r4V, ...
    x>2+r2w/2, r5V); %fplot(y)
a = linspace(-5,5,1000);
b = double(subs(y,a)); % Convert the symbolic function y(f) a vector

plot(a,b)
grid on

% eof