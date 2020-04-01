%% Written by Steven Bopp January 14, 2016
%% Angle Between Two Planes
% Changing the values of the a, b and c tensors will give the angle between the two planes
a= [1 0 0]; % p-vector
b= [1 1 1]; % q-vector
m= [1 0 0; 0 1 0; 0 0 1]; % Cubic Metric Tensor

c=norm(a*m*b'); d=sqrt(norm(a*m*a')); e=sqrt(norm(b*m*b')); 
angle=acosd(c/(d*e)) % Final calculation

%% Elastic Modulus in a Certain Direction for Cubic Structures
E100=125; E111=273; % Known moduli
a= [1 1 0]; % Direction in which you wnat to find a modulus
x= [1 0 0]; y= [0 1 0]; z= [0 0 1]; % All independant directions
m= [1 0 0; 0 1 0; 0 0 1]; % Cubic Metric Tensor

alpha= acosd(norm(a*m*x')/(norm(a*m*a')*norm(x*m*x'))); % Metric tensor angle calculations
beta=  acosd(norm(a*m*y')/(norm(a*m*a')*norm(y*m*y')));
gamma= acosd(norm(a*m*z')/(norm(a*m*a')*norm(z*m*z')));

%Final calculation
Euvw= (((1/E100)-3*(1/E100 - 1/E111)) * (alpha^2*beta^2 + beta^2*gamma^2 + gamma^2*alpha^2))^(-1)