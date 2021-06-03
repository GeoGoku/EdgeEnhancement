% INSTRUCTIONS FOR RUNNING NHF
% 
% NHF.m is a Matlab program for performing the edge enhancement described in the accompanying paper. 
% Place the program in your work directory along with the demonstration data file
% Pass the MatLab 2020a test

clc;
clear all;

% 'parallelepiped.mat' coressponds to section 'The single parallelepiped'
% 'magData.mat' and 'noisyMagData.mat' coressponds to section 'Three prisms'
% 'bishop.mat' coressponds to section 'The Bishop model'
obs = load('magData.mat');
obs = obs.obs;
x = obs.x;
y = obs.y;
v = obs.v;

figure;
pcolor(y, x, v); shading interp; hold on;
xlabel('East (m)');
ylabel('North (m)');
title('potential field data')
colorbar;

% now we can get the response of NHF (normalized Harris filter)
f = 1e-8;   % adjust f to get desired result.
k = 1;   % we recommend not to change this
[NHFR, R, upperEnv, ind1, ind2] = NHF(x, y, v, f, k);

figure;
pcolor(y, x, NHFR); shading interp; hold on;
xlabel('East (m)');
ylabel('North (m)');
title('NHF response')
colorbar;

% E-mail me if you have any problems.
% Tao Chen
% chentaosx@hotmail.com